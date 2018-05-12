#ifndef __UTL_H_INCLUDED
#define __UTL_H_INCLUDED

#define UTL_USE_ONLY_CPP_HEADERS



/*===========================================================================*/
/* Linear Algebra utilities
/*---------------------------------------------------------------------------*/

// This section may be redefined in other projects.
#if !defined(__CPPUTILS_H_INCLUDED__) && !defined(UTL_USE_ONLY_CPP_HEADERS)
#include <Eigen/Dense>
#include <math.h>
#include <vector>
#include <fstream>
namespace la
{

  template<typename T,unsigned int N>
  struct vec_t
  {
    typedef Eigen::Matrix<T,N,1> type;
  };

  typedef vec_t<double,2>::type dvec2_t;
  typedef vec_t<double,3>::type dvec3_t;
  typedef vec_t<double,4>::type dvec4_t;

  typedef vec_t<float,2>::type fvec2_t;
  typedef vec_t<float,3>::type fvec3_t;
  typedef vec_t<float,4>::type fvec4_t;

  typedef vec_t<int,2>::type ivec2_t;
  typedef vec_t<int,3>::type ivec3_t;
  typedef vec_t<int,4>::type ivec4_t;

  typedef vec_t<unsigned int,2>::type uivec2_t;
  typedef vec_t<unsigned int,3>::type uivec3_t;
  typedef vec_t<unsigned int,4>::type uivec4_t;

  template<typename T,unsigned int N>
  inline double tri_area(const typename vec_t<T,N>::type &p,
                    const typename vec_t<T,N>::type &q,
                    const typename vec_t<T,N>::type &r)
  {
    double a = (p-q).norm();
    double b = (q-r).norm();
    double c = (r-p).norm();

    double s = (a+b+c)/2.0;

    double area2  = s*(s-a)*(s-b)*(s-c);

    double area =  sqrt(area2);

    return area;
  }

  template<typename T,unsigned int R, unsigned int C>
  struct mat_t
  {
    typedef Eigen::Matrix<T,R,C> type;
  };

  typedef mat_t<double,3,3>::type dmat3x3_t;
  typedef mat_t<double,4,3>::type dmat4x3_t;
  typedef mat_t<double,4,4>::type dmat4x4_t;


}
#endif
/*---------------------------------------------------------------------------*/


#if !defined(UTL_USE_ONLY_CPP_HEADERS)
namespace utl {
/*---------------------------------------------------------------------------*/

/// \brief Read a matrix from binary file
///
/// \note No endian conversion is done

template<typename _S, int _R, int _C, int _O, int _MR, int _MC> inline
bool readBinMat(std::istream &file, Eigen::Matrix<_S,_R,_C,_O,_MR,_MC> & mat)
{
  // Data files are in ColMajor, read into a colmajor copy and copy it to mat
  if(_O & Eigen::RowMajor)
  {
    Eigen::Matrix<_S, _R, _C, Eigen::ColMajor, _MR, _MC> matIn;

    bool ret = readBinMat(file,matIn);
    if(ret) mat =matIn;
    return ret;
  }


  long int rc[2];
  file.read((char*)(void*)rc,sizeof(long int)*2);
  mat.resize(rc[0],rc[1]);
  file.read((char*)(void*)&mat(0,0),rc[0]*rc[1]*sizeof(_S));
  return true;
}

/*---------------------------------------------------------------------------*/

/// \brief Write a matrix to binary file
///
/// \note No endian conversion is done

template<typename _S, int _R, int _C, int _O, int _MR, int _MC> inline
bool writeBinMat(std::ostream &file, const Eigen::Matrix<_S,_R,_C,_O,_MR,_MC> & mat)
{
  // We only write in Col major, so we write using a colMajor copy
  if(_O & Eigen::RowMajor)
  {
    Eigen::Matrix<_S, _R, _C, Eigen::ColMajor, _MR, _MC> matOut=mat;
    return writeBinMat(file,matOut);
  }

  long int rc[] = {mat.rows(),mat.cols()} ;
  file.write((char*)(void*)rc,sizeof(long int)*2);
  file.write((char*)(void*)&mat(0,0),rc[0]*rc[1]*sizeof(_S));

  return true;
}

/*---------------------------------------------------------------------------*/

template<typename _S, int _R, int _C, int _O, int _MR, int _MC> inline
bool readBinMat(std::string file, Eigen::Matrix<_S, _R, _C, _O, _MR, _MC> &mat)
{
  std::fstream fs(file.c_str(),std::ios::binary|std::ios::in);
  if(!fs.is_open()) return false;
  return readBinMat(fs,mat);
}

/*---------------------------------------------------------------------------*/

template<typename _S, int _R, int _C, int _O, int _MR, int _MC> inline
bool writeBinMat(std::string file, const Eigen::Matrix<_S, _R, _C, _O, _MR, _MC> &mat)
{
  std::fstream fs(file.c_str(),std::ios::binary|std::ios::out);
  if(!fs.is_open()) return false;
  return writeBinMat(fs,mat);
}
/*---------------------------------------------------------------------------*/

/// \brief Base64 decode
void decode64(char * dat, int nbytes, const std::string &val);

/*---------------------------------------------------------------------------*/

/// \brief Base64 encode
std::string encode64(const char * dat, int n);

/*---------------------------------------------------------------------------*/

/// \brief Convert a typename to string
template <typename T> struct TypeToName {
  /// \brief Converter function
  static const char* Get(){return typeid(T).name();}
};
/*---------------------------------------------------------------------------*/
}
/*===========================================================================*/
#endif



/*===========================================================================*/
/* Misc utility functions
/*---------------------------------------------------------------------------*/
#include <tuple>
#include <iostream>
#include <vector>
namespace std {
	namespace detail {

		template <size_t n, typename... T>
		typename std::enable_if<(n >= sizeof...(T))>::type
			print_tuple(std::ostream&, const std::tuple<T...>&)
		{}

		template <size_t n, typename... T>
		typename std::enable_if<(n < sizeof...(T))>::type
			print_tuple(std::ostream& os, const std::tuple<T...>& tup)
		{
			if (n != 0)
				os << ", ";
			os << std::get<n>(tup);
			print_tuple<n + 1>(os, tup);
		}

		template <size_t n, typename... T>
		typename std::enable_if<(n >= sizeof...(T))>::type
			read_tuple(std::istream&, std::tuple<T...>&)
		{}

		template <size_t n, typename... T>
		typename std::enable_if<(n < sizeof...(T))>::type
			read_tuple(std::istream& is, std::tuple<T...>& tup)
		{
			is >> std::get<n>(tup);
			read_tuple<n + 1>(is, tup);
		}
	}

	template <typename... T>
	std::istream& operator>>(std::istream& is, std::tuple<T...>& tup)
	{
		detail::read_tuple<0>(is, tup);	return is;
	}

	template <typename... T>
	std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
	{
		detail::print_tuple<0>(os, tup); return os;
	}
}

namespace utl {

/*---------------------------------------------------------------------------*/

///  \brief A generic to string impl that uses the << operator of type T
template <typename T>
inline std::string to_string (const T& v)
{std::stringstream ss; ss << v ; return ss.str();}

/*---------------------------------------------------------------------------*/

template <class iter_t>
void argsort(iter_t b, iter_t e, std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

template <class iter_t>
void rearrange(iter_t b, iter_t e,const std::vector<size_t>& idxs);

/*---------------------------------------------------------------------------*/

/// \brief strip front and back whitespaces in given string
std::string trim(const std::string &s);

/*---------------------------------------------------------------------------*/

/// \brief Converts a given string to the template type
template <class T> inline T from_string(std::string s) {
	T v; std::stringstream(s) >> v; return v;
}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1>
inline void from_string(std::string s,T1& t1)
{std::stringstream ss(s);ss >> t1;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2>
inline void from_string(std::string s,T1& t1,T2& t2)
{std::stringstream ss(s);ss >> t1 >> t2;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2,typename T3>
inline void from_string(std::string s,T1& t1,T2& t2,T3& t3)
{std::stringstream ss(s);ss >> t1 >> t2 >> t3;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2,typename T3,typename T4>
inline void from_string(std::string s,T1& t1,T2& t2,T3& t3,T4& t4)
{std::stringstream ss(s);ss >> t1 >> t2 >> t3 >> t4;}

/*---------------------------------------------------------------------------*/

///  \brief Converts a given string to the template types
template <typename T1,typename T2,typename T3,typename T4,typename T5>
inline void from_string(std::string s,T1& t1,T2& t2,T3& t3,T4& t4, T5& t5)
{std::stringstream ss(s);ss >> t1 >> t2 >> t3 >> t4 >> t5;}


/*---------------------------------------------------------------------------*/

///  \brief String begins with
inline bool beginswith(std::string s_str, std::string c_str)
{return (s_str.substr(0,c_str.size()) == c_str);}

/*---------------------------------------------------------------------------*/

///  \brief String ends with
inline bool endswith(std::string s, std::string c)
{return (s.substr(s.size() - std::min(c.size(),s.size()),c.size()) == c);}

/*---------------------------------------------------------------------------*/

/// \brief Split a string into a vector of strings
void split(const std::string &s, const char* delim, std::vector<std::string> & v);

/*---------------------------------------------------------------------------*/

/// \brief Split a string into a vector of strings
inline std::vector<std::string> split(const std::string &s, const char* delim)
{
	std::vector<std::string> v; split(s, delim, v); return v;
}

}// namespace utl
/*===========================================================================*/






/*===========================================================================*/
/* Misc utility classes
/*---------------------------------------------------------------------------*/


#include <memory>
#include <iostream>
#include <thread>
#include <chrono>
#include <mutex>
#include <ctime>
#include <atomic>


#if !defined(UTL_USE_ONLY_CPP_HEADERS)
#include <boost/iterator/iterator_facade.hpp>
namespace utl {

/*---------------------------------------------------------------------------*/


/**
	\brief A class to iterate over non-comment lines of a file

	\note All parts of a line that succeed the comment char ('#' default) are
	assumed to be a comment and are stripped

	\note Sample usage
	for(file_line_iterator lgen("file.txt"),lend; lgen != lend; )
	cout << *lgen++;

**/

class file_line_iterator
	: public boost::iterator_facade<
	file_line_iterator
	, std::string const
	, boost::forward_traversal_tag
	>
{
public:
	file_line_iterator(const char * f, char c_char = '#');
	file_line_iterator(std::string f, char c_char = '#') :file_line_iterator(f.c_str(), c_char) {}
	file_line_iterator() {}

private:
	friend class boost::iterator_core_access;
	void increment();
	bool equal(file_line_iterator const& other) const;
	const std::string &dereference() const;

private:
	std::shared_ptr<std::ifstream> is;
	std::string                      value;
	char                             c_char;
};
}

#endif

namespace utl {

/*---------------------------------------------------------------------------*/

///** \brief Simple stopwatch timer (to measure wall clock time NOT cpu time) **/
//class timer
//{
// public:
//  timer()
//  {restart();}
//
//  inline void   restart()
//  { _start_time = boost::posix_time::microsec_clock::local_time(); }
//
//  inline double elapsed() const
//  {
//    boost::posix_time::time_duration td =
//        boost::posix_time::microsec_clock::local_time() - _start_time;
//
//    return double(td.total_milliseconds())/1000;
//  }
//
// private:
//  boost::posix_time::ptime _start_time;
//}; // timer

/*---------------------------------------------------------------------------*/

/** \brief A simple multi-threaded logger                                   **/
class logger
{
 public:

  enum severity_level {trace,debug,info,warning,error,fatal};

  inline void push(const std::string & log)
  {

	std::time_t cur_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	
	std::lock_guard<std::mutex>  lock(s_mutex);

	std::clog << " [" << std::string(std::ctime(&cur_time)).substr(0, 19) << "]"
		<< std::hex
		<< " [ 0x" << std::this_thread::get_id() << "]"
		<< std::dec << log << std::endl;
  }

  static inline logger& get() {return s_logger;}

private:

  static std::mutex   s_mutex;
  static logger       s_logger;
}; // logger

namespace detail
{
/** \brief a small extension of std::pair where only first is initialized   **/
template<class _T1, class _T2> struct pair:public std::pair<_T1,_T2>
{pair(const _T1 &t1){std::pair<_T1,_T2>::first = t1;}};
}//namespace detail


// \brief a simple c++11 atomics based mutex 
class atomic_mutex_t {
	std::atomic_flag flag = ATOMIC_FLAG_INIT;
public:

	inline void lock() { while (flag.test_and_set(std::memory_order_acquire)); }
	inline void unlock() { flag.clear(std::memory_order_release); }
};

// \brief a basic concurrent vector implementation. 
// Its a bit incomplete for now and methods that are not reimplemented here
// Might cause concurrency issues. Use with caution in mt environments
template<typename T> class concurrent_vector :protected std::vector<T> {

public:

	typedef std::vector<T> base_t;
	typedef atomic_mutex_t mutex_t;

	inline void push_back(const typename base_t::value_type& _Val) {
		std::lock_guard<mutex_t> l(mutex);
		base_t::push_back(_Val);
	}

	inline void clear() {
		std::lock_guard<mutex_t> l(mutex);
		base_t::clear();
	}

	inline void reserve(size_t n) {
		std::lock_guard<mutex_t> l(mutex);
		base_t::reserve(n);
	}

	inline void resize(size_t n) {
		std::lock_guard<mutex_t> l(mutex);
		base_t::resize(n);
	}

protected:

	mutex_t mutex;

};

}// namespace utl

/*---------------------------------------------------------------------------*/
#if !defined(UTL_USE_ONLY_CPP_HEADERS)
namespace utl {

	// \brief Unordered pair type, for use as keys in maps/ sets etc
	// EG : utl::uoPair_t<int> a(1,0), b(0,1), c(1,1);
	//      if(a == b ) std::cout  << "a == b" << std::endl; 
	//      if(a != c ) std::cout  << "a != c" << std::endl; 
	//
	// Will result in output as :
	// a == b
	// a != b

	template <typename T>
	struct uoPair_t : std::pair<T, T> {

		typedef std::pair<T, T> base_t;
		typedef uoPair_t         this_t;


		constexpr uoPair_t() : base_t() {}
		constexpr uoPair_t(const T& _Val1, const T& _Val2) : base_t(_Val1, _Val2) {}
		constexpr uoPair_t(const uoPair_t&) = default;
		constexpr uoPair_t(uoPair_t&&) = default;
		constexpr uoPair_t(const base_t&a) :base_t(a) {};
		constexpr uoPair_t(base_t&& a) :base_t(a) {};
#ifdef MSVC
		template<class _Other1, class _Other2,
		class = typename enable_if<is_convertible<const _Other1&, T>::value
			&& is_convertible<const _Other2&, T>::value, void>::type>
			constexpr uoPair_t(const pair<_Other1, _Other2>& _Right)
			: first(_Right.first), second(_Right.second)
		{}
#else
                template<class _U1, class _U2, class = typename
                         std::enable_if<std::__and_<std::is_convertible<const _U1&, T>,
                                          std::is_convertible<const _U2&, T>>::value>::type>
                  constexpr uoPair_t(const std::pair<_U1, _U2>& __p)
                  : base_t(__p.first,__p.second){ }

#endif
		template<class _Other1, class _Other2>
		this_t& operator=(const std::pair<_Other1, _Other2>& _Right) {
			base_t::operator=(_Right);return (*this);
		}

		friend bool operator==(const uoPair_t &a, const uoPair_t &b) {
			bool ret = ((a.first == b.first && a.second == b.second) || (a.first == b.second && b.first == a.second));
			// std::cout <<"operator==( [" << a.first <<","  << a.second <<"],[" << b.first <<","  << b.second <<"])=" << ret<< std::endl;		
			return ret;
		}
	};

	typedef uoPair_t<int> uoIntPair_t;
	typedef uoPair_t<std::string> uoStrPair_t;

	inline std::string getDirName(std::string path) {
		size_t found = path.find_last_of("/\\");

		if (found == std::string::npos) return "";
		else return (path.substr(0, found) + "/");
	}

	inline std::string getFileName(std::string path) {
		size_t found = path.find_last_of("/\\");

		if (found == std::string::npos) return path;
		else return path.substr(found+1);
	}


}

namespace std {
	template<typename T>
	struct hash<utl::uoPair_t<T> > {	// hash functor for uoPair_t

		size_t operator()(const utl::uoPair_t<T>& x) const
		{
			size_t hash = (x.first < x.second) ?
				(std::hash<T>()(x.first) ^ std::hash<T>()(x.second)) :
				(std::hash<T>()(x.second) ^ std::hash<T>()(x.first));
			// std::cout <<"Hash(" << x.first <<","  << x.second <<")=" << hash<< std::endl;
			return hash;
		}
	};
}

#endif

namespace std {
	// Clamp will be available in C++ 17

	template<class T, class Compare>
	constexpr const T& clamp(const T& v, const T& lo, const T& hi, Compare comp)
	{
#ifdef MSVC
		return assert(!comp(hi, lo)),
			comp(v, lo) ? lo : comp(hi, v) ? hi : v;
#else
	  return comp(v,lo) ? lo : comp(hi, v) ? hi : v;
#endif
	}

	template<class T>
	constexpr const T& clamp(const T& v, const T& lo, const T& hi)
	{
		return clamp<T>(v, lo, hi, std::less<T>());
	}

}

/*===========================================================================*/




/*===========================================================================*/
/* Macro definitions
/*---------------------------------------------------------------------------*/

#define is_in_range(i,b,e) (((b) <= (i)) && ((i) < (e)))

/*---------------------------------------------------------------------------*/


#ifdef WIN32
#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#else
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define LOG_MACRO(__isEnabled)  \
  if(__isEnabled)\
  for(utl::detail::pair<bool,std::stringstream> __utl_lm_v__(true); \
      __utl_lm_v__.first ;__utl_lm_v__.first=false,\
  utl::logger::get().push(__utl_lm_v__.second.str())) \
	__utl_lm_v__.second << "["<<__FILENAME__<<","<<__func__<<","<<__LINE__<<"]"



#ifdef UTL_ENABLE_TRACELOGS
#define TLOG LOG_MACRO(true)
#else
#define TLOG LOG_MACRO(false)
#endif

#ifndef NDEBUG
#define DLOG LOG_MACRO(true)
#else
#define DLOG LOG_MACRO(false)
#endif
#define ILOG LOG_MACRO(true)
#define WLOG LOG_MACRO(true)
#define ELOG LOG_MACRO(true)
#define FLOG LOG_MACRO(true)

/*---------------------------------------------------------------------------*/

#define SVAR(VAR) " "<<#VAR << " = "<< (VAR)

/*---------------------------------------------------------------------------*/

#ifndef ASSERT
#ifndef NDEBUG
#define ASSERT(cond) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERT(cond)
#endif // ifndef NDEBUG
#endif // ifndef ASSERT

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV
#ifndef NDEBUG
#define ASSERTV(cond,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV(cond,var1)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV2
#ifndef NDEBUG
#define ASSERTV2(cond,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV2(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#ifndef ASSERTV3
#ifndef NDEBUG
#define ASSERTV3(cond,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to assert condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n"\
  ;throw std::runtime_error(ss.str());}
#else  //ifndef NDEBUG
#define ASSERTV3(cond,var1,var2)
#endif // ifndef NDEBUG
#endif // ifndef ASSERTV2

/*---------------------------------------------------------------------------*/

#define ENSURE(cond,mesg) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV(cond,mesg,var1) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV2(cond,mesg,var1,var2) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#define ENSUREV3(cond,mesg,var1,var2,var3) if (!(cond))\
{ std::stringstream ss; \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "\
    <<"Message : " << #mesg << "\n"\
    << #var1 <<" = "<< (var1) << "\n"\
    << #var2 <<" = "<< (var2) << "\n"\
    << #var3 <<" = "<< (var3) << "\n";\
  ;throw std::runtime_error(ss.str());}

/*---------------------------------------------------------------------------*/

#if defined(WIN32)
void __dump_error__(const std::string & s);
#else
inline void __dump_error__(const std::string & s) 
{std::cerr << s << std::endl;throw std::runtime_error(s);}
#endif

#define ENSURES(cond) \
  if(!(cond)) \
  for(std::stringstream ss ; true ; __dump_error__(ss.str())) \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "

#define THISCODEPOINT \
	(std::string("(") + (__FILE__) + ","+__func__+","+std::to_string(__LINE__)+")")

/*---------------------------------------------------------------------------*/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/*---------------------------------------------------------------------------*/
#define BIT(x) (1<<(x))
/*---------------------------------------------------------------------------*/
#define FWDDECLARE_STRUCT(TYPE) \
struct TYPE; \
typedef std::shared_ptr<TYPE> TYPE ## Ptr;\
typedef std::shared_ptr<const TYPE> TYPE ## CPtr; \


/*---------------------------------------------------------------------------*/

namespace utl {
	namespace __tscout_detail {

		struct tsguard {

			tsguard() :guard(tscout_mutex) {}

			static utl::atomic_mutex_t tscout_mutex;

			std::lock_guard<utl::atomic_mutex_t> guard;

			bool done = false;

		};
	}
}

// Threadsafe cout
#define TSCOUT for(utl::__tscout_detail::tsguard i; !i.done; i.done = true) \
std::cout<< "0x"<<std::hex<<std::this_thread::get_id()<<": "

//#define tscout for(utl::__tscout_detail::tsguard i; !i.done; i.done = true) \
//std::cout<< "0x"<<std::hex<<std::this_thread::get_id()<<": "




/*===========================================================================*/


#endif
