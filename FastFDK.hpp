#ifndef FASTFDK_H_
#define FASTFDK_H_

/*
 * Author : Christophe Meessen
 *          meessen@cppm.in2p3.fr
 *
 *   Date : 2 jun 2009
 */
 

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "ims_float.hpp"

#ifndef DEBUG
#undef assert
#define assert(A)
#endif

using namespace std;

// define computation precision
typedef double Real;

class FastFDK
{
public:

	/*! @brief Construct FastFDK object
	 *  
	 * @param iw image width in pixels
	 * @param ih image height in pixels
     * @param id number of images in scan
	 * @param ps pixel size in mm
	 * @param yd distance from origin to detector in mm
	 * @param sd distance from source to origin in mm 
	 * @param vw volume width in voxels
	 * @param vd volume depth in voxels
	 * @param vs voxel size in mm
	 * @param na angular span of scan in degree
	 */
	FastFDK( uint16_t const iw, uint16_t const ih, uint16_t const id, 
				Real const ps, Real const yd, Real const sd, uint16_t const vw, 
				uint16_t const vd, Real const vs, Real const na );

	//! cleanup internal storage
	~FastFDK();
	
	//! reconstruct volume from images
	ims_float &reconstruct( const ims_float &images );

	//! return pixel value in image at coordinates u and v or 0 if out of bounds
	float getPixel( float const *img, int16_t const u, int16_t const v );
	
	//! filter image
	void filter( float *outImg, const float* inImg );

	/*! @brief backproject image into volume at angle a
	 *
	 * @param img array of pixel values of image
	 * @param a angle of volume for retroprojection
	 */
	void backproject( float const *img, Real const a, uint16_t const n );

	//! normalize voxel values
	void normalize( float epsilon = 1e-4 );

	/*! @brief bilinear interpolation of values using dx and dy as ratio
	 *  
	 * @param dx distance ratio between v1 and v2 or v3 and v4
	 * @param dy distance ratio between v1 and v3 or v2 and v4
	 * @param v1 top left value
	 * @param v2 top right value
	 * @param v3 bottom left value
	 * @param v4 bottom right value
	 */
	static float binLinInterp( Real const dx, Real const dy, 
				float const v1, float const v2, float const v3, float const v4 );

	//! return reference on volume
	ims_float &getVolume() { return m_vol; }

protected:
    const int16_t m_iw; ///< image width in pixels
    const int16_t m_ih; ///< image height in pixels
    const int16_t m_id; ///< number of images of scan
    Real m_ps;          ///< pixel size in mm
    const Real m_yd;    ///< distance from origin to detector
    const Real m_sd;    ///< distance from source to origin
    const uint16_t m_vw; ///< volum width in voxels
    const uint16_t m_vd; ///< volum depth in voxels
    const Real m_vs;    ///< voxel size in mm
    const Real m_na;    ///< angular span of scan in degree
    float *m_img;       ///< storage for filtered image
    bool  *m_cyl;       ///< flag set to true if voxel in cylinder
    ims_float m_vol;    ///< output volume
    float *m_w0;        ///< initial image weightening
    float *m_buf;       ///< storage for image filtering
    float *m_flt;       ///< odd values of convolution filter
    Real m_k;           ///< flt[0]
    uint16_t m_nf;      ///< number of filter values

};

// ----------------- INLINE -----------------

// constructor initializes precomputed data
inline FastFDK::FastFDK( uint16_t const iw, uint16_t const ih, uint16_t const id, 
				Real const ps, Real const yd, Real const sd, uint16_t const vw, 
				uint16_t const vd, Real const vs, Real const na ) 
	: m_iw(iw), m_ih(ih), m_id(id), m_ps(ps), m_yd(yd), m_sd(sd), 
		m_vw(vw), m_vd(vd), m_vs(vs), m_na(na), m_vol( vw, vw, vd ),
		m_w0(NULL), m_buf(NULL), m_flt(NULL)
{
	// assume images contain the rotation axis and origin

	// adjust pixel size so that distance from source to detector is now m_sd
	m_ps *= m_sd/(m_sd+m_yd);  
	
	// initialize initial image weighting
	m_w0 = new float[m_iw*m_ih];

	// for all pixels i,j compute u,v coordinates and pixel weighting
	float *p = m_w0;
	Real const k0 = m_sd*m_sd;
	for( int j = 0; j < m_ih; ++j )
	{
		Real const v = (m_ih/2-0.5-j)*m_ps;
		Real const k1 = k0 + v*v;
		for( int i = 0; i < m_iw; ++i )
		{
			Real u = (i+0.5-m_iw/2)*m_ps;
			*p++ = float(m_sd/::sqrt( u*u + k1 ));
		}
	}

	// instantiate buffer storage for convolution
	m_buf = new float[m_iw];
	
	// instantiate buffer for filtered images
	m_img = new float[iw*ih*id];

	// initialize filtre values (generate only half of it)
	m_nf = m_iw/2;
	m_flt = new float[m_nf];
	memset( m_flt, 0, m_nf*sizeof(float) );

	// flt[0] filter value
	m_k = 1/(4*m_ps*m_ps);

	// compute odd values of filter
	Real const k2 = (Real)(M_PI*m_ps);
	Real sum = 0;
	for( uint16_t i = 0, j = 1; i < m_nf; ++i, j += 2 )
	{
		Real const tmp = k2*j;
		sum += m_flt[i] = -1/(tmp*tmp);
	}
	//m_k = -2*sum;
	
	// precompute 'in cylinder' flag
	m_cyl = new bool[m_vw*m_vw];
	{
		bool *p = m_cyl;
		for( int i = 0; i < m_vw; ++i )
			for( int j = 0; j < m_vw; ++j, ++p )
				*p = ((2*i-m_vw)*(2*i-m_vw) + (2*j-m_vw)*(2*j-m_vw)) < (m_vw*m_vw);  
	}
}

// cleanup internal storage
inline FastFDK::~FastFDK()
{
	delete[] m_w0;
	delete[] m_buf;
	delete[] m_flt;
	delete[] m_img;
	delete[] m_cyl;
}

//! reconstruct volume from images
inline ims_float &FastFDK::reconstruct( const ims_float &images )
{
	// clear volume
	m_vol.fill( 0 );

	// filter all images
	uint32_t imgSize = m_iw*m_ih;
	float const *inImg = images.pixels();
	float *outImg = m_img;
	for( int i = 0; i < m_id; ++i, inImg += imgSize, outImg += imgSize )
		filter( outImg, inImg );

	// compute angle step in radian (rotate in reverse direction)
	Real const angStep = -m_na*M_PI/(m_id*180);
	
	// sequentially back project filtered images
	Real ang = 0;
	inImg = m_img;
	for( int i = 0; i < m_id; ++i, inImg += imgSize, ang += angStep )
		backproject( inImg, ang, i+1 );

	normalize();
	return m_vol;
}

// return pixel value in image at coordinates u and v or 0 if out of bounds
inline float FastFDK::getPixel( float const *img, int16_t const u, int16_t const v )
{
	return (u >= 0 && u < m_iw && v >= 0 && v < m_ih )? img[u + m_iw*v] : 0;
}

// filter image
inline void FastFDK::filter( float *outImg, const float* inImg )
{
	float const *wp = m_w0, *ip = inImg;
	float *op = outImg;
	for( uint16_t j = 0; j < m_ih; ++j, ip += m_iw, wp += m_iw, op += m_iw )
	{
		// copy image row into buffer and apply weighting filter
		memcpy( m_buf, ip, m_iw*sizeof(float) );
		for( uint16_t i = 0; i < m_iw; ++i )
			m_buf[i] *= wp[i];

		// apply convolution filter
		for( uint16_t i = 0; i < m_iw; ++i )
		{
			Real v0 = m_buf[i]*m_k;
			for( uint16_t k = 0, l = 1; k < m_nf; ++k, l += 2 )
			{
				if( (i - l) >= 0 )
				{
					if( (i + l) < m_iw )
						v0 += (m_buf[i-l] + m_buf[i+l]) * m_flt[k];
					else
						v0 += m_buf[i-l] * m_flt[k];
				}
				else if( (i + l) < m_iw )
					v0 += m_buf[i+l] * m_flt[k];
				else
					break;
			}
			op[i] = (float)v0;
		}
	}
}

// back project image of specified angle into volume
inline void FastFDK::backproject( float const *img, Real const a, 
	uint16_t const n )
{
	// precompute trigonometric value
	const Real ca = ::cos( a ), sa = ::sin( a );
	
	// precompute detector projection ratio
	const Real k1 = m_sd/m_ps;
	
	// precompute image size and volume slice size
    //uint32_t const imgSize = m_iw*m_ih;
	uint32_t const volSize = m_vw*m_vw;

	// for horizontal voxels rows 
	for( uint32_t j = 0, l = 0; j < m_vw; ++j )
	{
		// compute y coordinate component of voxel center
		Real const vy = (m_vw/2 - 0.5 - j)*m_vs;
		
		// for all voxels in row
		for( uint16_t i = 0; i < m_vw; ++i, ++l )
		{
			// skip voxels out of cylinder
			if( !m_cyl[l] )
				continue;
				
			// compute x coordinate component of voxel center
			Real const vx = (i + 0.5 - m_vw/2)*m_vs;

			// apply volume rotation
			Real const vvx = vx * ca - vy * sa;
			Real const vvy = vy * ca + vx * sa;

			// compute projection scaling factor
			Real const k2 = k1/(vvy+m_sd);

			// compute u coordinates on detector in pixel units: (0,0) is image center
			Real const u = vvx*k2 + m_iw/2 - 0.5;
			
			// compute right pixel column index (use nearest integer rounding)
			int16_t iu = (uint16_t)u;
			if( u < 0 )
				--iu;						
			
			// compute remainder
			Real const du = u - iu;

			// sanity check
			assert( du >= 0 && du <= 1 );
			
			// for all voxels in the volume column
			float *pvxl = m_vol.pixels() + l;
			for( int k = 0; k < m_vd; ++k, pvxl += volSize )
			{
				// compute z coordinate component of voxel center
				Real const vz = (m_vd/2 - 0.5 - k)*m_vs;

				// compute v coordinates on detector in pixel units
				Real const v = m_ih/2 - vz*k2 - 0.5;
				
				// compute top pixel row index (use nearest integer rounding)
				int16_t iv = (uint16_t)v;
				if( v < 0 )
					--iv;
				
				// compute remainder so that v = iv + dv
				Real const dv = v - iv;

			/*
				cout << " vx: " << vvx << " vy: " << vvy << " vz: " << vz; 
				cout << " u: " << u << " iu: " << iu << " du: " << du;
				cout << " v: " << v << " iv: " << iv << " dv: " << dv;
				cout << endl;
			*/
			
				// sanity check
				assert( dv >= 0 && dv <= 1 );
			
				// get pixel value using bilinear interpolation
				double val = binLinInterp( du, dv, 
					getPixel( img, iu, iv ), getPixel( img, iu+1, iv ), 
					getPixel( img, iu, iv+1 ), getPixel( img, iu+1, iv+1 ) ); 

				// compute post projection weightening
				//Real k3 = m_sd/(m_sd+vy*ca-vx*sa);
				Real k3 = m_sd/(m_sd+vvy);
				val *= k3 * k3; 
									// accumulate averaged filtered back projected image into volume
				*pvxl += ((float)val - *pvxl)/n;
			}
		}
	}
}

// normalize voxel values
inline void FastFDK::normalize( float epsilon )
{
	for( float *p = m_vol.pixels(), *e = p + m_vol.stack_size(); p != e; ++p )
		if( *p < epsilon && *p != 0 )
			*p = 0;
		else
			//*p /= 2*M_PI; // ?
			*p /= 2; // ?
}

// bilinear interpolation of values using dx and dy as ratio
inline float FastFDK::binLinInterp( Real const dx, Real const dy, 
		float const v1, float const v2, float const v3, float const v4 )
{
	assert( dx >= 0 && dx <= 1 );
	assert( dy >= 0 && dy <= 1 );
	Real const k1 = v1+dx*(v2-v1), k2 = v3+dx*(v4-v3);
	return float(k1 + dy*(k2 - k1)); 
}

#endif

