#ifndef ACTS_VALUE_CORRECTOR_H
#define ACTS_VALUE_CORRECTOR_H 1

// ACTS include(s)
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts
{
  /// @cond detail
  namespace detail
  {
    /**
     * @brief check and correct parameter values
     *
     * Values in the given vector are interpreted as values for the given parameters. As those
     * they are checked whether they are inside the allowed range and corrected if necessary.
     *
     * Invocation:
     *   - `value_corrector<params...>::result(parVector)` where `parVector` contains
     *     `sizeof...(params)` elements
     *
     * @post All values in the argument `parVector` are within the valid parameter range.
     *
     * @tparam params template parameter pack containing the multiple parameter identifiers
     */
    template<ParID_t... params>
    struct value_corrector;

    /// @cond
    template<typename R,ParID_t... params>
    struct value_corrector_impl;

    template<ParID_t... params>
    struct value_corrector
    {
      typedef ActsVector<ParValue_t,sizeof...(params)> ParVector_t;

      static void result(ParVector_t& values)
      {
        value_corrector_impl<ParVector_t,params...>::calculate(values,0);
      }
    };

    template<typename R, ParID_t first,ParID_t... others>
    struct value_corrector_impl<R,first,others...>
    {
      static void calculate(R& values,unsigned int pos)
      {
        typedef typename par_type<first>::type parameter_type;
        if(parameter_type::may_modify_value)
          values(pos) = parameter_type::getValue(values(pos));
        value_corrector_impl<R,others...>::calculate(values,pos+1);
      }
    };

    template<typename R,ParID_t last>
    struct value_corrector_impl<R,last>
    {
      static void calculate(R& values,unsigned int pos)
      {
        typedef typename par_type<last>::type parameter_type;
        if(parameter_type::may_modify_value)
          values(pos) = parameter_type::getValue(values(pos));
      }
    };
    /// @endcond
  }  // end of namespace detail
  /// @endcond
}  // end of namespace Acts

#endif // ACTS_VALUE_CORRECTOR_H