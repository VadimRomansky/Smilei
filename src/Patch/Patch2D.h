#ifndef PATCH2D_H
#define PATCH2D_H

#include "Patch.h"
#include "Field2D.h"

class SimWindow;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch2D : public Patch
{
public:
    //! Constructor for Patch
    Patch2D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved);

    //! Destructor for Patch
    virtual ~Patch2D() {};


    // MPI exchange/sum methods for particles/fields
    //   - fields communication specified per geometry (pure virtual)
    // --------------------------------------------------------------

    //! init comm / sum densities
    virtual void initSumField( Field* field, int iDim );
    //! finalize comm / sum densities
    virtual void finalizeSumField( Field* field, int iDim );

    //! init comm / exchange fields
    virtual void initExchange( Field* field );
    //! finalize comm / exchange fields
    virtual void finalizeExchange( Field* field );
    //! init comm / exchange fields in direction iDim only
    virtual void initExchange( Field* field, int iDim );
    //! finalize comm / exchange fields in direction iDim only
    virtual void finalizeExchange( Field* field, int iDim );

    // Create MPI_Datatype to exchange fields
    virtual void createType( Params& params );

    //! MPI_Datatype to sum [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    //! MPI_Datatype to exchange [ndims_+1][iDim=0 prim/dial][iDim=1 prim/dial]
    //!   - +1 : an additional type to exchange clrw lines
    MPI_Datatype ntype_[3][2][2];
    // Use a buffer per direction to exchange data before summing
    Field2D buf[2][2];

    //! Should be pure virtual, see Patch class
    inline bool isWestern()  { return locateOnBorders(0, 0); }
    //! Should be pure virtual, see Patch class
    inline bool isEastern()  { return locateOnBorders(0, 1); }
    //! Should be pure virtual, see Patch class
    inline bool isSouthern() { return locateOnBorders(1, 0); }
    //! Should be pure virtual, see Patch class
    inline bool isNorthern() { return locateOnBorders(1, 1); }

    //! Return MPI rank of this->hrank +/- 1
    //! Should be replaced by an analytic formula
    inline int getMPIRank(int hrank_pm1) {
	if  (hrank_pm1 == neighbor_[0][0]) return MPI_neighbor_[0][0];
	else if  (hrank_pm1 == neighbor_[0][1]) return MPI_neighbor_[0][1];
	else if  (hrank_pm1 == neighbor_[1][0]) return MPI_neighbor_[1][0];
	else if  (hrank_pm1 == neighbor_[1][1]) return MPI_neighbor_[1][1];
	else
	    return MPI_PROC_NULL;
    }



};

#endif