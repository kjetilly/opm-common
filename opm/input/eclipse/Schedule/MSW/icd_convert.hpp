#ifndef OPM_ICD_CONVERT_HPP
#define OPM_ICD_CONVERT_HPP

namespace Opm {

template<typename T>
T from_int(long long int_status);

template<typename T>
long long to_int(T status);

}

#endif
