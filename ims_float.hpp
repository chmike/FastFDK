#ifndef _IMS_FLOAT_H_
#define _IMS_FLOAT_H_

/* Author: Christophe Meessen <meessen@cppm.in2p3.fr>
   Revision: 2009/02/17
   Version: 0.3
*/

/*! @brief Image Stack class with pixels of type float 

   This class is derived from the ims_base class and provide support
   for adding a segment holding some basic stats information on the
   pixel values. 
 
   The image stack type is identifid by the value 1 and the pixel size
   is 4 bytes. 

   The segment holding the stats is identified by the ASCII letter 
   sequence 'F', 'l', 't' and 'S' and hold four values. 
   The minimum and maximum float values in the stack stored as floats.
   The mean and variance value stored as double precison floats.

   The method update_stats() will compute the statistic values and
   add the required segment if not present. If pixel valus are 
   modified, it is required to call update_stats() to recompute the
   statistic values.

   The load and save method will call update_stats() themselves if
   the segment is not present which is the case for new image stacks. 
*/

#include "ims_base.hpp"

#include <iostream>

// define float ims base class, add statistic information
class ims_float : public ims_base<float,0,sizeof(float)>
{
protected:
	//! define base type 
	typedef ims_base<float,0,sizeof(float)> base_type; 

public:
	//! declare type name for segmentMap
	typedef base_type::segmentMap segmentMap;

	//! default constructor
	ims_float() : m_min(0), m_max(0), m_mean(0), m_var(0) {}

	/*! @brief instantiating a new initialized ims of specified size
	 *
	 *  @param iw image width (pixel)
	 *  @param ih image height(pixel)
	 *  @param id stack depth (images)
	 *  @param initVal initial pixel value
	 */
	ims_float( int iw, int ih, int id, const float initVal = 0 ) : base_type( iw, ih, id, initVal ),
		m_min(initVal), m_max(initVal), m_mean(initVal), m_var(0) {}

	/*! @brief copy initializer
	 *
	 *  @param ims Image stack to copy
	 */
	ims_float( const ims_float &ims ) : base_type( ims )
	{
		// copy stats
		m_min = ims.m_min;
		m_max = ims.m_max;
		m_mean = ims.m_mean;
		m_var = ims.m_var;
	}

	/*! @brief construct with memory loading specified file
	 *
	 *  @param file name to load as memory mapped file stack
	 */
	ims_float( const char* fileName ) : base_type( fileName )
	{
		// copy float stats in member variables if present
		load_stats();
	}

	/*! @brief copy ims into current instance overriding its content
	 *  @param ims image stack to copy
	 */
	void assign( const ims_float &ims )
	{
		// call parent assignement method
		base_type::assign( ims );

		// copy stats
		m_min = ims.m_min;
		m_max = ims.m_max;
		m_mean = ims.m_mean;
		m_var = ims.m_var;
	}

	/*! @brief assignment operator	
	 *
	 *  @param ims_float Image stack to assign
	 */
	ims_float& operator=( const ims_float &ims ) { assign( ims ); return *this; }

	//! recompute stats and add FltS segment if required
	void update_stats()
	{
		// compute stats over pixels in a single pass
		float *data = pixels();
		m_max = m_min = *data;
		m_mean = m_var = 0;
		uint64_t n = 0;
		for( uint64_t i = 0; i < m_sSz; ++i )
		{
			float val = data[i];
			if( val < m_min )
				m_min = val;
			else if( val > m_max )
				m_max = val;
			++n;
			double delta = val - m_mean;
			m_mean += delta/n;
			m_var += delta*( val - m_mean );
		}
		m_var /= n - 1;

		// locate float stats segment in image stack
		const uint32_t FltS = getFltS();
		segmentMap::iterator it = m_seg.find( FltS );

		// if it is missing create the segment
		if( it == m_seg.end() )
		{
			uint32_t sz = sizeof(float)*2 + sizeof(double)*2;
			m_seg.insert( std::make_pair( FltS, segment( FltS, sz, new uint8_t[sz], true ) ) );
			it = m_seg.find( FltS );
			m_hdr.sn++; 
		}

		// copy member variable data in segment data
		float *dataFlt = (float*)it->second.data;
		dataFlt[0] = m_min;
		dataFlt[1] = m_max;
		double *dataDbl = (double*)(dataFlt+2);
		dataDbl[0] = m_mean;
		dataDbl[1] = m_var;
	}

	//! return minimum value
	float min_val() const { return m_min; }

	//! return maximum value
	float max_val() const { return m_max; }

	//! return mean value
	double mean() const { return m_mean; }

	//! return variance value
	double variance() const { return m_var; }

	/*! @brief Saves an ims_float file
	 *			
	 *  throws if failed to save the ims file.
	 *
	 *  @param fileName file to save as ims file
	 */    
	void save( const char* fileName ) const
	{
		// make sure we have float stats in image
		if( m_seg.find( getFltS() ) == m_seg.end() )
			const_cast<ims_float*>(this)->update_stats();

		// call parent save method
		base_type::save( fileName );
	}

	/*! @brief Load image stack into memory
	 *
	 * Pixel data may be modified. Won't modify original file. 
	 *
	 *  @param fileName file to open with memory mapping
	 */
	void load( const char* fileName )
	{
		// call parent load method
		base_type::load( fileName );

		// copy float stats in member variables if present
		load_stats();
	}

	//! @brief return float stat segment id
	unsigned long getFltS() const { return 0x53746C46UL; }

protected:

	//! copy stats value from segment to member variables
	void load_stats()
	{
		// locate float stats segment in image stack
		segmentMap::iterator it = m_seg.find( getFltS() );

		// if found then copy that into member variables
		if( it != m_seg.end() )
		{
			float *dataFlt = (float*)it->second.data;
			m_min = dataFlt[0];
			m_max = dataFlt[1];
			double *dataDbl = (double*)(dataFlt+2);
			m_mean = dataDbl[0];
			m_var = dataDbl[1];
		}
		else
			update_stats();
	}

	float m_min;   //< minimum float value
	float m_max;   //< maximum float value
	double m_mean; //< mean value
	double m_var;  //< variance
};

#endif
