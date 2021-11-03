#ifndef DIAGNOSTICPARTICLEBINNINGBASE_H
#define DIAGNOSTICPARTICLEBINNINGBASE_H

#include "Diagnostic.h"

#include "Histogram.h"

class DiagnosticParticleBinningBase : public Diagnostic
{
    friend class SmileiMPI;
    friend class VectorPatch;
    
public :

    //! Default constructor
    DiagnosticParticleBinningBase(
        Params &params,
        SmileiMPI *smpi,
        Patch *patch,
        int diagId,
        std::string diagName,
        bool time_accumulate,
        PyObject *deposited_quantity,
        std::vector<std::string> excluded_Axes
    );
    //! Default destructor
    ~DiagnosticParticleBinningBase() override;
    
    virtual void openFile( Params &params, SmileiMPI *smpi ) override;
    
    void closeFile() override;
    
    bool prepare( int itime ) override;
    
    void calculate_auto_limits( Patch *patch, SimWindow *simWindow, unsigned int ipatch );
    
    virtual void run( Patch *patch, int itime, SimWindow *simWindow ) override;
    
    virtual bool writeNow( int itime );
    
    void write( int itime, SmileiMPI *smpi ) override;
    
    //! Clear the array
    virtual void clear();
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override
    {
        int size = output_size*sizeof( double );
        // + data_array + index_array +  axis_array
        // + nparts_max * (sizeof(double)+sizeof(int)+sizeof(double))
        return size;
    };
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint( int istart, int istop, Patch *patch ) override;
    
    std::vector<std::vector<double> > patches_mins, patches_maxs;
    
protected:
    
    //! True for Screen only
    bool time_accumulate;
    
    //! number of timesteps during which outputs are averaged
    int time_average;
    
    //! list of the species that will be accounted for
    std::vector<unsigned int> species_indices;
    
    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;
    
    //! Histogram object
    Histogram *histogram;
    
    unsigned int output_size;
    
    int total_axes;
    std::vector<hsize_t> dims;
    
//    //! Minimum and maximum spatial coordinates that are useful for this diag
//    std::vector<double> spatial_min, spatial_max;
};

#endif

