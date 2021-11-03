#include "PartWall.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "Species.h"
#include "BoundaryConditionType.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;

// Partwall constructor
PartWall::PartWall( double pos, unsigned short dir, string kind, double dt ) :
    position( pos ),
    direction( dir )
{
    // Define the "wall" function pointer
    if( kind == "reflective" ) {
        wall = &reflect_particle_wall;
    } else if( kind == "remove" ) {
        wall = &remove_particle_wall;
    } else if( kind == "stop" ) {
        wall = &stop_particle_wall;
    } else if( kind == "thermalize" ) {
        wall = &thermalize_particle_wall;
    }

    dt_ = dt;
   
}

// Applies the wall's boundary condition to one particle
void PartWall::apply( Species *species, int imin, int imax, std::vector<double> &invgf, Random * rand, double &energy_change )
{
    ( *wall )( species, imin, imax, direction, position, dt_, invgf, rand, energy_change );
}


// Reads the input file and creates the ParWall objects accordingly
PartWalls::PartWalls( Params &params, Patch *patch )
{
    
    resize( 0 );
    unsigned int numpartwall=PyTools::nComponents( "PartWall" );
    if( patch->isMaster() && numpartwall>0 ) {
        MESSAGE( 1, "Adding particle walls:" );
    }
    direction.resize( numpartwall );
    position .resize( numpartwall );
    kind     .resize( numpartwall );
    
    // Loop over each wall component and parse info
    for( unsigned int iwall = 0; iwall < numpartwall; iwall++ ) {
    
        // Extract the direction of the wall
        string dirstring;
        bool extract_x = PyTools::extractOrNone( "x", position[iwall], "PartWall", iwall );
        bool extract_y = PyTools::extractOrNone( "y", position[iwall], "PartWall", iwall );
        bool extract_z = PyTools::extractOrNone( "z", position[iwall], "PartWall", iwall );
        if( extract_x + extract_y + extract_z > 1 ) {
            ERROR( "PartWall #" << iwall << ": cannot have several locations (x, y or z)" );
        }
        if( extract_x ) {
            direction[iwall] = 0;
            dirstring = "x";
        } else if( extract_y ) {
            if( params.nDim_particle < 2 ) {
                ERROR( "PartWall #" << iwall << ": cannot have y-location in 1D" );
            }
            direction[iwall] = 1;
            dirstring = "y";
        } else if( extract_z ) {
            if( params.nDim_particle < 3 ) {
                ERROR( "PartWall #" << iwall << ": cannot have z-location y in 1D or 2D" );
            }
            direction[iwall] = 2;
            dirstring = "z";
        } else {
            ERROR( "PartWall #" << iwall << " must have one location (x, y or z)" );
        }
        
        // Ewtract the kind of wall
        PyTools::extract( "kind", kind[iwall], "PartWall", iwall );
        if( kind[iwall].empty() || ( kind[iwall]!="reflective" && kind[iwall]!="remove" && kind[iwall]!="stop" && kind[iwall]!="thermalize" ) ) {
            ERROR( "For PartWall #" << iwall << ", `kind` must be one of reflective, remove, stop, thermalize" );
        }
        
        // Find out wether this proc has the wall or not
        if( position[iwall] >= patch->getDomainLocalMin( direction[iwall] )
                && position[iwall] <= patch->getDomainLocalMax( direction[iwall] ) ) {
            push_back( new PartWall( position[iwall], direction[iwall], kind[iwall], params.timestep ) );
        }
        
        // Create new wall
        MESSAGE( 2, "Adding a wall at "<<dirstring<<" = "<< position[iwall] << ", kind:" << kind[iwall] << ( kind[iwall]=="thermalize" ? " thermCond" : "" ) );
    }
}

// Clones an existing vector of partWalls
PartWalls::PartWalls( PartWalls *partWalls, Patch *patch )
{
    // Copy all the walls info, so that all patches know about all walls
    resize( 0 );
    direction = partWalls->direction;
    position  = partWalls->position ;
    kind      = partWalls->kind     ;
    
    // Create walls, but only those within the current domain
    unsigned int nwalls=direction.size();
    for( unsigned int iwall = 0; iwall < nwalls; iwall++ ) {
        if( position[iwall] >= patch->getDomainLocalMin( direction[iwall] )
                && position[iwall] <= patch->getDomainLocalMax( direction[iwall] ) ) {
            push_back( new PartWall( position[iwall], direction[iwall], kind[iwall], partWalls->vecPartWall[iwall]->dt_ ) );
        }
    }

}


// Destructor
PartWalls::~PartWalls()
{
    int nwalls=size();
    for( int i=0; i<nwalls; i++ ) {
        delete vecPartWall[i];
    }
    clear();
}
