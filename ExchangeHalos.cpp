// Example Includes
#include "ExchangeHalos.hpp"

// C++ includes
#include <iostream>
#include <string>

static double evaluate_function( double lon, double lat, int type = 1, double multiplier = 1.0 )
{
    switch( type )
    {
        case 1:
            return ( 2.0 + std::pow( sin( 2.0 * lat ), 16.0 ) * cos( 16.0 * lon ) ) * multiplier;
        default:
            return ( 2.0 + cos( lon ) * cos( lon ) * cos( 2.0 * lat ) ) * multiplier;
    }
}

moab::ErrorCode RuntimeContext::create_sv_tags( moab::Tag& tagScalar, moab::Tag& tagVector,
                                                moab::Range& entities ) const
{
    // Get element (centroid) coordinates so that we can evaluate some arbitrary data
    std::vector< double > entCoords = compute_centroids( entities );  // [entities * [lon, lat]]

    if( proc_id == 0 ) std::cout << "> Getting scalar tag handle " << scalar_tagname << "..." << std::endl;
    double defSTagValue = -1.0;
    bool createdTScalar = false;
    // Create the scalar exchange tag: default name = "scalar_variable"
    runchk( moab_interface->tag_get_handle( scalar_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, tagScalar,
                                            moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, &defSTagValue, &createdTScalar ),
            "Retrieving scalar tag handle failed" );

    assert( createdTScalar );
    // set the data for scalar tag
    {
        std::vector< double > tagValues( entities.size(), -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = index++ * 2;
            return evaluate_function( entCoords[offset], entCoords[offset + 1] );
        } );
        // Set local scalar tag data for exchange
        runchk( moab_interface->tag_set_data( tagScalar, entities, tagValues.data() ),
                "Setting scalar tag data failed" );
    }

    if( proc_id == 0 ) std::cout << "> Getting vector tag handle " << vector_tagname << "..." << std::endl;
    std::vector< double > defVTagValue( vector_length, -1.0 );
    bool createdTVector = false;
    // Create the scalar exchange tag: default name = "vector_variable"
    runchk( moab_interface->tag_get_handle( vector_tagname.c_str(), vector_length, moab::MB_TYPE_DOUBLE, tagVector,
                                            moab::MB_TAG_CREAT | moab::MB_TAG_DENSE, defVTagValue.data(),
                                            &createdTVector ),
            "Retrieving vector tag handle failed" );

    assert( createdTVector );
    // set the data for vector tag
    {
        const int veclength = vector_length;
        std::vector< double > tagValues( entities.size() * veclength, -1.0 );
        std::generate( tagValues.begin(), tagValues.end(), [=, &entCoords]() {
            static int index = 0;
            const int offset = ( index++ / veclength ) * 2;
            return evaluate_function( entCoords[offset], entCoords[offset + 1], 2, ( index % veclength + 1.0 ) );
        } );
        // Set local tag data for exchange
        runchk( moab_interface->tag_set_data( tagVector, entities, tagValues.data() ),
                "Setting vector tag data failed" );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode RuntimeContext::load_file( bool load_ghosts )
{
    /// Parallel Read options:
    ///   PARALLEL = type {READ_PART} : Read on all tasks
    ///   PARTITION_METHOD = RCBZOLTAN : Use Zoltan partitioner to compute an online partition and redistribute on the fly
    ///   PARTITION = PARALLEL_PARTITION : Partition as you read based on part information stored in h5m file
    ///   PARALLEL_RESOLVE_SHARED_ENTS : Communicate to all processors to get the shared adjacencies
    ///   consistently in parallel
    ///   PARALLEL_GHOSTS : a.b.c
    ///                   : a = 2 - highest dimension of entities (2D in this case)
    ///                   : b = 1 - dimension of entities to calculate adjacencies (vertex=0, edges=1)
    ///                   : c = 3 - number of ghost layers needed (3 in this case)
    std::string read_options   = "DEBUG_IO=0;";
    std::string::size_type idx = input_filename.rfind( '.' );
    std::string extension      = "";
    if( num_procs > 1 && idx != std::string::npos )
    {
        extension = input_filename.substr( idx + 1 );
        if( !extension.compare( "nc" ) )
            // PARTITION_METHOD= [RCBZOLTAN, TRIVIAL]
            read_options += "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;"
                            "PARALLEL_RESOLVE_SHARED_ENTS;NO_EDGES;NO_MIXED_ELEMENTS;VARIABLE=;";
        else if( !extension.compare( "h5m" ) )
            read_options +=
                "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                "PARALLEL_RESOLVE_SHARED_ENTS;" +
                ( load_ghosts ? "PARALLEL_THIN_GHOST_LAYER;PARALLEL_GHOSTS=2.1." + std::to_string( ghost_layers ) + ";"
                              : "" );
        else
        {
            std::cout << "Error unsupported file type (only h5m and nc) for this example: " << input_filename
                      << std::endl;
            return moab::MB_UNSUPPORTED_OPERATION;
        }
    }

    // Load the file from disk with given read options in parallel
    return moab_interface->load_file( input_filename.c_str(), &fileset, read_options.c_str() );
}

std::vector< double > RuntimeContext::compute_centroids( const moab::Range& entities ) const
{
    double node[3];
    std::vector< double > eCentroids( entities.size() * 2 );  // [lon, lat]
    size_t offset = 0;
    for( auto entity : entities )
    {
        // Get the element coordinates (centroid) on the real mesh
        runchk_cont( moab_interface->get_coords( &entity, 1, node ), "Getting entity coordinates failed" );

        // scale by magnitude so that mesh is on unit sphere
        double magnitude = std::sqrt( node[0] * node[0] + node[1] * node[1] + node[2] * node[2] );
        node[0] /= magnitude;
        node[1] /= magnitude;
        node[2] /= magnitude;

        // compute the spherical transformation onto unit sphere
        eCentroids[offset] = atan2( node[1], node[0] );
        if( eCentroids[offset] < 0.0 ) eCentroids[offset] += 2.0 * M_PI;
        eCentroids[offset + 1] = asin( node[2] );

        offset += 2;  // increment the offset
    }
    return eCentroids;
}
