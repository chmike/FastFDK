#ifndef _IMS_BASE_H_
#define _IMS_BASE_H_

/* Author: Christophe Meessen <meessen@cppm.in2p3.fr>
   Revision: 2009/02/17
   Version: 0.4
*/

/*! @brief Image Stack template class 

    This class may be used to create, load or save image stacks with
    associated information called segments. An image stack is a 
    sequence of images or a 3D volume as produced by tomodensitometry.

    Image stacks are saved in files of type ".ims" with a format
    described below. The class provides a methode to load a stack of
    images in ims format and save a stack of images in this format.   

    The pixel type is defined as a template argument so that the same
    class can be used for different pixel types. The constrain is
    however that the pixel storage size is a multiple of a byte.
    In addition to the pixel class type, the integer value of its
    identification in ims files must be given with the byte size of
    a pixel. 

    IMS file format
    ---------------

    Graphical description.  

    |<------ 32 bits ------>|
    |_______________________|
    |     |     |     |     |
    | 'I' | 'M' | 'S' |  1  |
    |_____|_____|_____|_____|
    |                       |
    | pixel type identifier |
    |_______________________|  File Header
    |           |           |
    |   width   |  height   |
    |___________|___________|
    |           |           |
    |   depth   |  nbrSeg   |
    |___________|___________|________________
    |                       |
    |                       |
    =     Pixel values      =   Pixel data
    |                       |
    |_______________________|________________
    |                       |
    |  Segment type ident.  |
    |_______________________| Segment 0 header
    |                       |
    |   Segment byte size   |
    |_______________________|________________
    |                       |
    |                       |
    =     Segment data      = Segment 0 payload
    |                       |
    |_______________________|________________
    |                       |
    |  Segment type ident.  |
    |_______________________| Segment 1 header
    |                       |
    |   Segment byte size   |
    |_______________________|________________
    |                       |
    |                       |
    =     Segment data      = Segment 1 payload
    |                       |
    |_______________________|________________
    |                       |
    |        . . .          |

  The first three bytes are the ASCII value of the corresponding letters
  identifying the file type. The fourth value is the encoding version. 
  It is currently 1. 

  The next 32 bit unsigned integer value is the pixel type identifier. 
  (i.e. 0 for floats). The next 16 bits unsigned integer value are 
  respectively the images width, the images heights, the stack depth and 
  the number of  segments.

  The header is followed by the pixel values stored images after images and
  for each image, rows by rows. The first pixel in an image is the top left
  most the last is the bottom right pixel. 

  Segments, stored after the pixels values, are stored in sequence.
  Each segment has an 8 Byte header and is following by 0 or m data bytes. 
  The first four bytes of the header is an array of ASCII values identifying
  the type of segment. For instance 'T', 'E', 'X' and 'T', in that order,
  identifies a section of type text encoded in utf8. The next 4 bytes of 
  the header encodes the number of bytes of data that follows. 
  It is encoded as a 32bit unsigned integer as all values stored in the
  data part of the segment.  
*/

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <cerrno>
#include <cassert>
#include <iostream>
#ifdef _WIN32
#include <windows.h>
typedef char int8_t;
typedef short int16_t;
typedef int int32_t;
typedef long long int64_t;
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdint.h>
#include <string.h>
#endif


template <class TPixel, unsigned long PixelTypeNo, int TPixelSize >
class ims_base
{
public:
	
	//! define own type
	typedef ims_base<TPixel, PixelTypeNo, TPixelSize> my_type; 

	//! segment information structure
	struct segment
	{
		//! default constructor
		segment() : type(0), size(0), data(NULL), need_delete(false) {}

		//! initialize segment value
		segment( uint32_t type, uint32_t size, uint8_t *data, bool need_delete ) 
			: type(type), size(size), data(data), need_delete(need_delete) {}
		
		uint32_t type;    //< segment type
		uint32_t size;    //< segment size
		uint8_t *data;    //< data size
		bool need_delete; //< true if data needs to be deleted
	};

	//! vector of segment information
	typedef std::map< uint32_t, segment > segmentMap;

	//! vector of images
	typedef std::vector<TPixel*> imageVct;

	//! @brief default constructor instantiating an uninitialized ims 
	ims_base() 
		: m_iSz(0), m_sSz(0), m_pxl(NULL), m_map(NULL), m_fSz(0)
#ifdef _WIN32
		, m_file(0), m_fmap(0)
#endif
		{ clear(); }

	/*! @brief instantiating a new initialized ims of specified size
	 *
	 *  @param iw image width (pixel)
	 *  @param ih image height(pixel)
	 *  @param id stack depth (images)
	 *  @param initVal initial pixel value
	 */
	ims_base( int iw, int ih, int id, const TPixel &initVal = TPixel() ) 
		: m_iSz(0), m_sSz(0), m_pxl(NULL), m_map(NULL), m_fSz(0)
#ifdef _WIN32
		, m_file(0), m_fmap(0)
#endif
		{ init( iw, ih, id, initVal ); }

	/*! @brief copy initializer
	 *
	 *  @param ims Image stack to copy
	 */
	 ims_base( const my_type &ims ) 
		: m_iSz(0), m_sSz(0), m_pxl(NULL), m_map(NULL), m_fSz(0)
#ifdef _WIN32
		, m_file(0), m_fmap(0)
#endif
		{ assign( ims ); }

	/*! @brief construct with memory loading specified file
	 *
	 *  @param file name to load as memory mapped file stack
	 */
	ims_base( const char* fileName )
		: m_iSz(0), m_sSz(0), m_pxl(NULL), m_map(NULL), m_fSz(0)
#ifdef _WIN32
		, m_file(0), m_fmap(0)
#endif
		{ load( fileName ); }

	//! @brief Image stack destructor
	~ims_base() { clear(); }

	//! @brief Release all space used by image stack. 
	void clear()
	{
		if( is_mapped() )
			unload();
		else
			delete[] m_pxl;
		m_img.clear();
		m_pxl = NULL;
		m_hdr = Hdr();
		m_iSz = 0; 
		m_sSz = 0;
		for( typename segmentMap::iterator it = m_seg.begin(); it != m_seg.end(); ++it )
			if( it->second.need_delete )
				delete[] it->second.data;
		m_seg.clear(); 
	}

	//! fill images with specified value
	void fill( TPixel val )
	{
		for( TPixel *p = pixels(), *e = p + stack_size(); p != e; ++p )
			*p = val; 
	}

	/*! @brief (re) initialize the image stack with the given size and value
	 *  @param iw image width (pixel)
	 *  @param ih image height(pixel)
	 *  @param id stack depth (images)
	 *  @param initVal initial pixel value
	 */
	void init( int iw, int ih, int id, const TPixel &initVal = TPixel() )
	{
		// check parameter validity
		if( iw < 0 || iw > 65536 || ih < 0 || ih > 65535 || id < 0 || id > 65535 )
			throw std::invalid_argument( "Invalid image stack dimension" );

		// clear instance
		clear();
	
		// set image stack dimension
		m_hdr.iw = iw;
		m_hdr.ih = ih;
		m_hdr.id = id;
		m_iSz = iw * ih;
		m_sSz = m_iSz * id;
		
		// instantiate and clear image stack
#ifdef _WIN32
		assert( m_sSz < (1<<30) );
		m_pxl = new TPixel[(unsigned long)m_sSz];
#else
		m_pxl = new TPixel[m_sSz];
#endif
		for( int i = 0; i < m_sSz; ++i )
			m_pxl[i] = initVal;

		// instantiate image vector
		m_img.reserve( m_hdr.id );
		for( int i = 0; i < m_hdr.id; ++i )
			m_img.push_back( m_pxl + i*m_iSz );
	}

	/*! @brief copy ims into current instance overriding its content
	 *  @param ims image stack to copy
	 */
	void assign( const my_type &ims )
	{
		// ignore if assign to self
		if( &ims == this )
			return;

		// clear current instance
		clear();

		// copy stack size information
		m_hdr = ims.m_hdr;
		m_iSz = ims.m_iSz;
		m_sSz = ims.m_sSz;

		// instantiate and copy pixel values
#ifdef _WIN32
		assert( m_sSz < (1<<30) );
		m_pxl = new TPixel[(unsigned long)m_sSz];
#else
		m_pxl = new TPixel[m_sSz];
#endif
        for( uint32_t i = 0; i < m_sSz; ++i )
			m_pxl[i] = ims.m_pxl[i];

		// instantiate image vector
		m_img.reserve( m_hdr.id );
		for( int i = 0; i < m_hdr.id; ++i )
			m_img.push_back( m_pxl + i*m_iSz );

		// copy segment map
		for( typename segmentMap::iterator it = m_seg.begin(); it != m_seg.end(); ++it )
		{
			// if segment data is not empty copy it
			uint8_t *data = NULL;
			if( it->second.size != 0 )
			{
				assert( it->second.data != NULL );
				data = new uint8_t[it->second.size];
				memcpy( data, it->second.data, it->second.size );
			} 

			// append segment copy to segment list
			m_seg.insert( make_pair( it->second.type, 
				segment( it->second.type, it->second.size, data, data != NULL ) ) );
		}
	}

	//! Return image width (Pixels)
	int width() const { return m_hdr.iw; }

	//! Return image height (Pixels)
	int height() const { return m_hdr.ih; }

	//! Return stack depth (Images)
	int depth() const { return m_hdr.id; }

	//! Return image size (Pixels)
	int image_size() const { return m_iSz; }

	//! Return stack size (Pixels)
	long long stack_size() const { return m_sSz; }

	//! Return reference on image vector
	imageVct &images() { return m_img; }

	//! Return reference on image vector
	const imageVct &images() const { return m_img; }

	//! Return pointer on first pixel
	TPixel *pixels() { return m_pxl; }

	//! Return pointer on first pixel
	const TPixel *pixels() const { return m_pxl; }

	//! return true if image stack is mapped into memory
	bool is_mapped() const { return m_map != NULL; }
	
	//! return segment tag corresponding to string
	static uint32_t getType( const char *type )
	{
		uint32_t v = 1;
		if( *(uint8_t*)(&v) == 1 )
			v = type[0]|(type[1]<<8)|(type[2]<<16)|(type[3]<<24);
		else	 
			memcpy( &v, type, 4 );
        return v;
	}

	//! return segment tag corresponding to string
	static uint32_t getType( const std::string &type )
		{ return getType( type.c_str() ); }

	/*!@brief get reference on segment vector
	 *
	 * Append, delete or modify segments, but manage associated memory.
	 * Allocated data with new[] or delete it with delete[].
	 * set the need_delete field accordingly.
	 *
	 *  @return reference on segment list
	 */
	segmentMap& segments() { return m_seg; }

	/*! @brief assignment operator	
	 *
	 *  @param ims Image stack to assign
	 */
	my_type& operator=( const my_type &ims ) { assign( ims ); return *this; }

	/*! @brief Saves an ims file
	 *			
	 *  throws if failed to save the ims file.
	 *
	 *  @param fileName file to save as ims file
	 */    
	void save( const char* fileName ) const
	{
		// open file for writing
		std::ofstream fic( fileName, std::ios::binary );
		if( !fic )
		{
			std::stringstream s;
			s << "Failed opening file '" << fileName << "' for saving";
			throw std::invalid_argument( s.str() );
		}
		
		// make sure all segments are referenced
		const_cast<Hdr*>(&m_hdr)->sn = m_seg.size();	

		// get copy of header
		uint8_t hdr[16];
		m_hdr.save( hdr );

		// write copy of header
		if( !fic.write( (char*)hdr, 16 ) )
		{
			std::stringstream s;
			s << "Failed writing header into file '" << fileName << "'";
			throw std::invalid_argument( s.str() );
		}

		// attempt to write pixels data
		uint64_t sz = m_sSz*TPixelSize;
		const uint32_t maxBuf = 1 << 30;
		uint8_t *p = (uint8_t*)m_pxl;
		while( sz > maxBuf )
		{
			if( !fic.write( (char*)p, maxBuf ) )
			{
				std::stringstream s;
				s << "Failed writing pixels into file '" << fileName << "'";
				throw std::invalid_argument( s.str() );
			}
			sz -= maxBuf;
			p += maxBuf;
		}
		if( !fic.write( (char*)p, (uint32_t)sz ) )
		{
			std::stringstream s;
			s << "Failed writing pixels into file '" << fileName << "'";
			throw std::invalid_argument( s.str() );
		}

		// for each segment
		for( typename segmentMap::const_iterator it = m_seg.begin(); it != m_seg.end(); ++it )
		{
			// set segment header
			const segment &seg(it->second);
			unsigned char segHdr[8];
			Hdr::saveUInt4( segHdr, seg.type );
			Hdr::saveUInt4( segHdr+4, seg.size );

			// attempt to write segment header
			if( !fic.write( (char*)segHdr, 8 ) )
			{
				std::stringstream s;
				s << "Failed writing segment header into file '" << fileName << "'";
				throw std::invalid_argument( s.str() );
			}

			// attempt to write segment data
			if( !fic.write( (char*)seg.data, seg.size ) )
			{
				std::stringstream s;
				s << "Failed writing segment data into file '" << fileName << "'";
				throw std::invalid_argument( s.str() );
			}
		}
	}

	/*! @brief use copy one write memory mapping to load image stack
	 *
	 * Pixel data may be modified. Won't modify original file. 
	 *
	 *  @param fileName file to open with memory mapping
	 */
	void load( const char* fileName )
	{
		uint8_t *data = NULL;
#ifdef _WIN32
		m_file = CreateFile( fileName, GENERIC_READ|GENERIC_WRITE, 0,
			NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL| 
			FILE_FLAG_WRITE_THROUGH|FILE_FLAG_NO_BUFFERING|
			FILE_FLAG_RANDOM_ACCESS, NULL );
		if( m_file == INVALID_HANDLE_VALUE )
		{
			m_file = 0;
			std::stringstream s;
			s << "Failed opening file '" << fileName << "': ";
			char buf[1024];
			strerror_s( buf, 1024, GetLastError() );
			s << buf;
			throw std::invalid_argument( s.str() );
		}

		m_fmap = CreateFileMapping( m_file, NULL, PAGE_READONLY, 0, 0, NULL );
		if( m_fmap == 0 )
		{
			CloseHandle( m_file );
			m_file = 0;
			std::stringstream s;
			s << "Failed mapping file '" << fileName << "': ";
			char buf[1024];
			strerror_s( buf, 1024, GetLastError() );
			s << buf;
			throw std::invalid_argument( s.str() );
		}
		
		data = (Hdr*)MapViewOfFile( m_fmap, FILE_MAP_COPY,0,0,0 );
		if( data == NULL )
		{ 
			CloseHandle( m_fmap );
			CloseHandle( m_file );
			m_fmap = m_file = 0;
			std::stringstream s;
			s << "Failed mapping file '" << fileName << "': ";
			char buf[1024];
			strerror_s( buf, 1024, GetLastError() );
			s << buf;
			throw std::invalid_argument( s.str() );
		}

		unsigned long low, high;
		low = GetFileSize( m_file, &high );
		m_fSz = high;
		m_fSz <<= 32;
		m_fSz |= low;
#else
		// try opening file
		//int fd = ::open( fileName, O_RDONLY|O_LARGEFILE );
		int fd = ::open( fileName, O_RDONLY );
		if( fd == -1 )
		{
			std::stringstream s;
			s << "Failed opening file '" << fileName << "': ";
			s << ::strerror( errno );
			throw std::invalid_argument( s.str() );
		}
 		
		// try getting file info and size
		struct stat st;
		if( fstat( fd, &st ) < 0 )
		{
			::close(fd);
			std::stringstream s;
			s << "Failed getting size of file '" << fileName << "': ";
			s << ::strerror( errno );
			throw std::invalid_argument( s.str() );
		}
		m_fSz = st.st_size;
		
		// map file in memory
		data = (uint8_t*)::mmap(NULL, st.st_size, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
		// close file any way
		::close(fd);

		// check header if mmap succeeded
		if( data == NULL )
		{
			m_fSz = 0;
			std::stringstream s;
			s << "Failed mapping file '" << fileName << "': ";
			s << ::strerror( errno );
			throw std::invalid_argument( s.str() );
		}

#endif
    // Instantiating header with header data
    Hdr hdr( data );

		// check for valid image stack check word and version
		if( hdr.chk != Hdr::CHECK_WORD )
		{
			m_map = data;
			unload();
			if( (hdr.chk ^ Hdr::CHECK_WORD)&0xFFFFFF )
			{
				std::stringstream s;
				s << "File '" << fileName << "' is not an image stack file";
				throw std::invalid_argument( s.str() );
			}
			std::stringstream s;
			s << "File '" << fileName << "' image stack version is ";
			s << (int)(hdr.chk >> 24) << " instead of ";
			s << (int)(Hdr::CHECK_WORD >> 24);
			throw std::invalid_argument( s.str() );
		}

		// clear current content
		clear();
		
		// initialize values
		m_map = data;
		m_hdr = hdr;
		m_iSz = m_hdr.iw*m_hdr.ih;
		m_sSz = m_iSz*m_hdr.id;
		m_pxl = (TPixel*)(data+16);

		// set image vector
		m_img.reserve( m_hdr.id );
		for( int i = 0; i < m_hdr.id; ++i )
			m_img.push_back( m_pxl + i*m_iSz );

		// setup segment map
		uint8_t *seg = (uint8_t*)(m_map + 16 + m_sSz*TPixelSize);
		for( int i = 0; i < m_hdr.sn; ++i )
		{
			// append segment copy to segment list
			uint32_t segTag = Hdr::loadUInt4( seg ), segSz = Hdr::loadUInt4( seg+4 );
			m_seg.insert( make_pair( segTag, segment( segTag, segSz, (uint8_t*)(seg+8), false ) ) );

			// jump to next segment
			seg += segSz + 8;
		}		
	}

	//! @brief return TEXT segment identifier
	unsigned long getTEXT() const { return 0x54584554UL; }

protected:

	//! @brief unloads the memory mapped file
	void unload()
	{
		// check that a file is mapped
		assert( is_mapped() );
#ifdef _WIN32
		UnmapViewOfFile( m_map );
		if( m_fmap )
			CloseHandle( m_fmap );
		if( m_file )
			CloseHandle( m_file );
		m_fmap = m_file = 0;
#else
		::munmap( m_map, m_fSz );
#endif
		m_map = NULL;
		m_fSz = 0;
	}

	//! Image Stack file header
	struct Hdr
	{
		//! default constructor
		Hdr() : chk(CHECK_WORD), it(PixelTypeNo), iw(0), ih(0), id(0), sn(0) {}

    //! construct with data stored at hdr
    Hdr( uint8_t *hdr ) { load( hdr ); }

		//! load header from buffer pointed by hdr
		void load( uint8_t *hdr )
		{
			chk = loadUInt4( hdr ); hdr += 4;
			it =  loadUInt4( hdr ); hdr += 4;
			iw =  loadUInt2( hdr ); hdr += 2;
			ih =  loadUInt2( hdr ); hdr += 2;
			id =  loadUInt2( hdr ); hdr += 2;
			sn =  loadUInt2( hdr );
		}
	
		//! write header content to cout
		void dump() const
		{
			std::cout << "chk: " << chk << std::endl;
			std::cout << "it : " <<  it << std::endl;
			std::cout << "iw : " <<  iw << std::endl;
			std::cout << "ih : " <<  ih << std::endl;
			std::cout << "id : " <<  id << std::endl;
			std::cout << "sn : " <<  sn << std::endl;
		}

		//! save header in buffer pointed by hdr
		void save( uint8_t *hdr ) const
		{
			saveUInt4( hdr, chk ); hdr += 4;
			saveUInt4( hdr, it  ); hdr += 4;
			saveUInt2( hdr, iw  ); hdr += 2;
			saveUInt2( hdr, ih  ); hdr += 2;
			saveUInt2( hdr, id  ); hdr += 2;
			saveUInt2( hdr, sn  );
		} 

		//! read two byte little endian unsigned integer stored at p
		static uint32_t loadUInt4( uint8_t *p )
			{ return p[0] | p[1]<<8 | p[2]<<16 | p[3]<<24; }

		//! read two byte little endian unsigned integer stored at p
		static uint16_t loadUInt2( uint8_t *p )
			{ return p[0] | p[1]<<8; }

		//! write four bytes little endian unsigned integer at p
		static void saveUInt4( uint8_t *p, uint32_t val )
		{
			*p++ = (uint8_t)val; val >>= 8; 
			*p++ = (uint8_t)val; val >>= 8; 
			*p++ = (uint8_t)val; val >>= 8; 
			*p++ = (uint8_t)val;
		}

		//! write two byte little endian unsigned integer at p
		static void saveUInt2( uint8_t *p, uint16_t val )
		{
			*p++ = (uint8_t)val; val >>= 8; 
			*p++ = (uint8_t)val;
		}

		//! "IMS/0"  little endian
		const static unsigned long CHECK_WORD = 0x01534D49; 
		unsigned long chk; //< check word and version field
		unsigned long  it; //< image type
		unsigned short iw; //< image width
		unsigned short ih; //< image height
		unsigned short id; //< image depth
		unsigned short sn; //< number of segments
  };

	Hdr         m_hdr;  //< image stack file header
	uint32_t    m_iSz;  //< image size   (pixel)  
	uint64_t    m_sSz;  //< stack size   (pixel)
	TPixel     *m_pxl;  //< Pointer on the first pixel
	imageVct    m_img;  //< vector of images
	uint8_t    *m_map;  //< pointer on mapped file (or NULL if unmapped)
	size_t      m_fSz;  //< mapped file size (0 if not mapped)
	segmentMap  m_seg;  //< segments mapped with type as key

#ifdef _WIN32
	HANDLE     m_file; //< handle on mapped file
	HANDLE     m_fmap; //< handle on file mapping
#endif
};


#endif
