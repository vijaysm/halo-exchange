#ifndef __ExchangeHalos_hpp_
#define __ExchangeHalos_hpp_

// MOAB includes
#include "moab/Core.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"

#ifndef MOAB_HAVE_MPI
#error "Please build MOAB with MPI..."
#endif

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

// C++ includes
#include <iostream>
#include <string>

#define dbgprint( MSG )                                           \
    do                                                            \
    {                                                             \
        if( context.proc_id == 0 ) std::cout << MSG << std::endl; \
    } while( false )

#define runchk( CODE, MSG )         \
    do                              \
    {                               \
        moab::ErrorCode err = CODE; \
        MB_CHK_SET_ERR( err, MSG ); \
    } while( false )

#define runchk_cont( CODE, MSG )                               \
    do                                                         \
    {                                                          \
        moab::ErrorCode err = CODE;                            \
        MB_CHK_ERR_CONT( err );                                \
        if( err ) std::cout << "Error:: " << MSG << std::endl; \
    } while( false )

/// @brief The RunttimeContext is an example specific class to store
/// the run specific input data, MOAB datastructures used during the run
/// and provides other utility functions to profile operations etc
struct RuntimeContext
{
  public:
    int dimension{ 2 };           /// dimension of the problem
    std::string input_filename;   /// input file name (nc format)
    std::string output_filename;  /// output file name (h5m format)
    int ghost_layers{ 3 };        /// number of ghost layers
    std::string scalar_tagname;   /// scalar tag name
    std::string vector_tagname;   /// vector tag name
    int vector_length{ 3 };       /// length of the vector tag components
    int num_max_exchange{ 10 };   /// total number of exchange iterations
    bool debug_output{ false };   /// write debug output information?
    int proc_id{ 1 };             /// process identifier
    int num_procs{ 1 };           /// total number of processes
    double last_counter{ 0.0 };   /// last time counter between push/pop timer

    // MOAB objects
    moab::Interface* moab_interface{ nullptr };
    moab::ParallelComm* parallel_communicator{ nullptr };
    moab::EntityHandle fileset{ 0 }, partnset{ 0 };

    /// @brief Constructor: allocate MOAB interface and communicator, and initialize
    /// other data members with some default values
    RuntimeContext( MPI_Comm comm = MPI_COMM_WORLD )
        : input_filename( std::string( MESH_DIR ) + std::string( "/io/mpasx1.642.t.2.nc" ) ),
          output_filename( "exchangeHalos_output.h5m" ), scalar_tagname( "scalar_variable" ),
          vector_tagname( "vector_variable" )
    {
        // Create the moab instance
        moab_interface = new( std::nothrow ) moab::Core;
        if( NULL == moab_interface ) exit( 1 );

        // Create sets for the mesh and partition.  Then pass these to the load_file functions to populate the mesh.
        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, fileset ), "Creating root set failed" );
        runchk_cont( moab_interface->create_meshset( moab::MESHSET_SET, partnset ), "Creating partition set failed" );

        // Create the parallel communicator object with the partition handle associated with MOAB
        parallel_communicator = moab::ParallelComm::get_pcomm( moab_interface, partnset, &comm );

        proc_id   = parallel_communicator->rank();
        num_procs = parallel_communicator->size();
    }

    /// @brief Destructor: deallocate MOAB interface and communicator
    ~RuntimeContext()
    {
        delete parallel_communicator;
        delete moab_interface;
    }

    /// @brief Parse the runtime command line options
    /// @param argc - number of command line arguments
    /// @param argv - command line arguments as string list
    void ParseCLOptions( int argc, char* argv[] )
    {
        ProgOptions opts;
        // Input mesh
        opts.addOpt< std::string >( "input", "Input mesh filename to load in parallel", &input_filename );
        // Output mesh
        opts.addOpt< void >( "debug", "Should we write output file? Default=false", &debug_output );
        opts.addOpt< std::string >(
            "output", "Output mesh filename for verification (use --debug). Default=exchangeHalos_output.h5m",
            &output_filename );
        // Dimension of the input mesh
        // Vector tag length
        opts.addOpt< int >( "vtaglength", "Size of vector components per each entity. Ddefault=3", &vector_length );
        // Number of halo (ghost) regions
        opts.addOpt< int >( "nghosts", "Number of ghost layers (halos) to exchange. Default=3", &ghost_layers );
        // Number of times to perform the halo exchange for timing
        opts.addOpt< int >( "nexchanges", "Number of ghost-halo exchange iterations to perform. Default=10",
                            &num_max_exchange );

        opts.parseCommandLine( argc, argv );
    }

    /// @brief Measure and start the timer to profile a task
    /// @param operation String name of the task being measured
    inline void timer_push( std::string operation )
    {
        mTimerOps = mTimer.time_since_birth();
        mOpName   = operation;
    }

    /// @brief Stop the timer and store the elapsed duration
    /// @param nruns Optional argument used to average the measured time
    void timer_pop( const int nruns = 1 )
    {
        double locElapsed = mTimer.time_since_birth() - mTimerOps;
        double avgElapsed = 0;
        double maxElapsed = 0;
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, parallel_communicator->comm() );
        MPI_Reduce( &locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, parallel_communicator->comm() );
        if( proc_id == 0 )
        {
            avgElapsed /= num_procs;
            if( nruns > 1 )
                std::cout << "[LOG] Time taken to " << mOpName.c_str() << ", averaged over " << nruns
                          << " runs : max = " << maxElapsed / nruns << ", avg = " << avgElapsed / nruns << "\n";
            else
                std::cout << "[LOG] Time taken to " << mOpName.c_str() << " : max = " << maxElapsed
                          << ", avg = " << avgElapsed << "\n";

            last_counter = maxElapsed / nruns;
        }
        mOpName.clear();
    }

    /// @brief Return the last elapsed time
    /// @return last_counter from timer_pop was called
    inline double last_elapsed() const
    {
        return last_counter;
    }

    /// @brief Load a MOAB supported file (h5m or nc format) from disk
    ///        representing an MPAS mesh
    /// @param load_ghosts Optional boolean to specify whether to load ghosts
    ///                    when reading the file (only relevant for h5m)
    /// @return Error code if any (else MB_SUCCESS)
    moab::ErrorCode load_file( bool load_ghosts = false );

    /// @brief Create scalar and vector tags in the MOAB mesh instance
    /// @param tagScalar Tag reference to the scalar field
    /// @param tagVector Tag reference to the vector field
    /// @param entities Entities on which both the scalar and vector fields are defined
    /// @return Error code if any (else MB_SUCCESS)
    moab::ErrorCode create_sv_tags( moab::Tag& tagScalar, moab::Tag& tagVector, moab::Range& entities ) const;

  private:
    /// @brief Compute the centroids of elements in 2D lat/lon space
    /// @param entities Entities to compute centroids
    /// @return Vector of centroids (as lat/lon)
    std::vector< double > compute_centroids( const moab::Range& entities ) const;

    moab::CpuTimer mTimer;
    double mTimerOps{ 0.0 };
    std::string mOpName;
};

#endif // #ifndef __ExchangeHalos_hpp_