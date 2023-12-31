/** @example ExchangeHalos Driver
 * \brief Example program that shows the use case for performing tag data exchange
 * between parallel processors in order to sync data on shared entities.
 *
 * <b>This example </b>:
 *    -# Initialize MPI and instantiates MOAB
 *    -# Gets user options: Input mesh file name, vector tag length, ghost layer size etc
 *    -# Create the root and partition sets
 *    -# Instantiate ParallelComm and read the mesh file in parallel using appropriate options
 *    -# Create the required number of ghost layers as requested by the user (default = 3)
 *    -# Get 2D MPAS polygonal entities in the mesh and filter to get only the "owned" entities
 *    -# Create two tags: scalar_variable (single data/cell) and vector_variable (multiple data/cell)
 *    -# Set tag data using analytical functions for both scalar and vector fields on owned entities
 *    -# Exchange shared entity information and tags between processors
 *      -# If debugging is turned on, store mesh file and tag on root process (will not contain data on shared entities)
 *      -# Perform exchange of scalar tag data on shared entities
 *      -# Perform exchange of vector tag data on shared entities
 *      -# If debugging is turned on, store mesh file and tag on root process (will now contain data on *all* entities)
 *    -#  Destroy the MOAB instance and finalize MPI
 *
 * <b>To run: </b>
 *      mpiexec -n np ./ExchangeHalos --input <mpas_mesh_file> --nghosts <ghostlayers> --vtaglength <vector component
 * size> \
 *                    --nexchanges <number of exchange runs>
 * <b>Example:</b>
 *      mpiexec -n 16 ./ExchangeHalos --input data/default_mesh_holes.h5m --nghosts 3 --vtaglength 100
 *
 * NOTE: --debug option can be added to write out extra files in h5m format to visualize some output (written from root
 * task only)
 *
 */
// Example Includes
#include "ExchangeHalos.hpp"

// C++ includes
#include <iostream>
#include <string>

using namespace moab;
using namespace std;

//
// Start of main test program
//
int main( int argc, char** argv )
{
    // Initialize MPI first
    MPI_Init( &argc, &argv );

    {
        // Create our context for this example run
        RuntimeContext context;
        dbgprint( "********** Exchange halos example **********\n" );

        // Get the input options
        context.ParseCLOptions( argc, argv );

        /////////////////////////////////////////////////////////////////////////
        // Print out the input parameters in use
        dbgprint( " -- Input Parameters -- " );
        dbgprint( "    Number of Processes  = " << context.num_procs );
        dbgprint( "    Input mesh           = " << context.input_filename );
        dbgprint( "    Ghost Layers         = " << context.ghost_layers );
        dbgprint( "    Scalar Tag name      = " << context.scalar_tagname );
        dbgprint( "    Vector Tag name      = " << context.vector_tagname );
        dbgprint( "    Vector Tag length    = " << context.vector_length << endl );
        /////////////////////////////////////////////////////////////////////////

        // Timer storage for all phases
        double elapsed_times[4];

        // Read the input file specified by user, in parallel, using appropriate options
        // Supports reading partitioned h5m files and MPAS nc files directly with online Zoltan partitioning
        context.timer_push( "Read input file" );
        {
            // Load the file from disk with given options
            runchk( context.load_file( false ), "MOAB::load_file failed for filename: " << context.input_filename );
        }
        context.timer_pop();
        elapsed_times[0] = context.last_elapsed();

        // Let the actual measurements begin...
        dbgprint( "\n- Starting execution -\n" );

        // We need to set up the ghost layers requested by the user. First correct for thin layers and then
        // call `exchange_ghost_cells` to prepare the mesh for use with halo regions
        context.timer_push( "Setup ghost layers" );
        {
            // Loop over the number of ghost layers needed and ask MOAB for layers 1 at a time
            for( int ighost = 0; ighost < context.ghost_layers; ++ighost )
            {
                // Exchange ghost cells
                int ghost_dimension  = context.dimension;
                int bridge_dimension = context.dimension - 1;
                // Let us now get all ghost layers from adjacent parts
                runchk( context.parallel_communicator->exchange_ghost_cells(
                            ghost_dimension, bridge_dimension, ( ighost + 1 ), 0, true /* store_remote_handles */,
                            true /* wait_all */, &context.fileset ),
                        "Exchange ghost cells failed" );  // true to store remote handles

                // Ensure that all processes understand about multi-shared vertices and entities
                // in case some adjacent parts are only m layers thick (where m < context.ghost_layers)
                if( ighost < context.ghost_layers - 1 )
                    runchk( context.parallel_communicator->correct_thin_ghost_layers(),
                            "Thin layer correction failed" );
            }
        }
        context.timer_pop();
        elapsed_times[1] = context.last_elapsed();

        // Get the 2D MPAS elements and filter it so that we have only owned elements
        Range dimEnts;
        {
            // Get all entities of dimension = dim
            runchk( context.moab_interface->get_entities_by_dimension( context.fileset, context.dimension, dimEnts ),
                    "Getting 2D entities failed" );
            // Get only owned entities! The ghosted/shared entities will get their data when we exchange
            // So let us filter entities based on the status: NOT x NOT_OWNED = OWNED status :-)
            runchk( context.parallel_communicator->filter_pstatus( dimEnts, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                    "Filtering pstatus failed" );

            // Aggregate the total number of elements in the mesh
            auto numEntities     = dimEnts.size();
            int numTotalEntities = 0;
            MPI_Reduce( &numEntities, &numTotalEntities, 1, MPI_INT, MPI_SUM, 0,
                        context.parallel_communicator->proc_config().proc_comm() );

            // We expect the total number of elements to be constant, immaterial of number of processes.
            // If not, we have a bug!
            dbgprint( "Total number of " << context.dimension << "D elements in the mesh = " << numTotalEntities );
        }

        Tag tagScalar = nullptr;
        Tag tagVector = nullptr;
        // Create two tag handles: scalar_variable and vector_variable
        // Set these tags with appropriate closed form functional data
        // based on element centroid information
        runchk( context.create_sv_tags( tagScalar, tagVector, dimEnts ), "Unable to create scalar and vector tags" );

        // let us write out the local mesh before tag_exchange is called
        // we expect to see data only on the owned entities - and ghosted entities should have default values
        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "> Writing to file *before* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_pre.h5m", "H5M", "" ),
                    "Writing to disk failed" );
        }

        // Perform exchange of tag data between neighboring tasks
        dbgprint( "> Exchanging tags between processors " );
        context.timer_push( "Exchange scalar tag data" );
        for( auto irun = 0; irun < context.num_max_exchange; ++irun )
        {
            // Exchange scalar tags between processors
            runchk( context.parallel_communicator->exchange_tags( tagScalar, dimEnts ),
                    "Exchanging scalar tag between processors failed" );
        }
        context.timer_pop( context.num_max_exchange );
        elapsed_times[2] = context.last_elapsed();

        context.timer_push( "Exchange vector tag data" );
        for( auto irun = 0; irun < context.num_max_exchange; ++irun )
        {
            // Exchange vector tags between processors
            runchk( context.parallel_communicator->exchange_tags( tagVector, dimEnts ),
                    "Exchanging vector tag between processors failed" );
        }
        context.timer_pop( context.num_max_exchange );
        elapsed_times[3] = context.last_elapsed();

        // let us write out the local mesh after tag_exchange is called
        // we expect to see real data on both owned and ghost entities in halo regions (non-default values)
        if( context.debug_output && ( context.proc_id == 0 ) )  // only on root process, for debugging
        {
            dbgprint( "> Writing to file *after* ghost exchange " );
            runchk( context.moab_interface->write_file( "exchangeHalos_output_rank0_post.h5m", "H5M", "" ),
                    "Writing to disk failed" );
        }

        // Write out the final mesh with the tag data and mesh -- just for verification
        if( context.debug_output )
        {
            dbgprint( "> Writing out the final mesh and data in MOAB h5m format. File = " << context.output_filename );
            string write_options = ( context.num_procs > 1 ? "PARALLEL=WRITE_PART;DEBUG_IO=0;" : "" );
            // Write out to output file to visualize reduction/exchange of tag data
            runchk( context.moab_interface->write_file( context.output_filename.c_str(), "H5M", write_options.c_str() ),
                    "File write failed" );
        }

        // Consolidated timing results: the data is listed as follows
        // [ntasks,  nghosts,  load_mesh(I/O),  exchange_ghost_cells(setup), exchange_tags(scalar),
        // exchange_tags(vector)]
        dbgprint( "\n> Consolidated: [" << context.num_procs << ", " << context.ghost_layers << ", " << elapsed_times[0]
                                        << ", " << elapsed_times[1] << ", " << elapsed_times[2] << ", "
                                        << elapsed_times[3] << "]," );

        // execution finished
        dbgprint( "\n********** ExchangeHalos Example DONE! **********" );
    }
    // Done, cleanup
    MPI_Finalize();

    return 0;
}
