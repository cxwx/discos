// -*- C++ -*-
//
// $Id: XBackendsErrorsC.h,v 1.1 2010/06/21 09:08:59 bliliana Exp $

// ****  Code generated by the The ACE ORB (TAO) IDL Compiler ****
// TAO and the TAO IDL Compiler have been developed by:
//       Center for Distributed Object Computing
//       Washington University
//       St. Louis, MO
//       USA
//       http://www.cs.wustl.edu/~schmidt/doc-center.html
// and
//       Distributed Object Computing Laboratory
//       University of California at Irvine
//       Irvine, CA
//       USA
//       http://doc.ece.uci.edu/
// and
//       Institute for Software Integrated Systems
//       Vanderbilt University
//       Nashville, TN
//       USA
//       http://www.isis.vanderbilt.edu/
//
// Information about TAO is available at:
//     http://www.cs.wustl.edu/~schmidt/TAO.html

// TAO_IDL - Generated from
// be/be_codegen.cpp:135

#ifndef _TAO_IDL____OBJECT_XBACKENDSERRORSC_H_
#define _TAO_IDL____OBJECT_XBACKENDSERRORSC_H_


#include /**/ "ace/config-all.h"

#if !defined (ACE_LACKS_PRAGMA_ONCE)
# pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */


#include "tao/AnyTypeCode/AnyTypeCode_methods.h"
#include "tao/ORB.h"
#include "tao/UserException.h"
#include "tao/Basic_Types.h"
#include "tao/String_Manager_T.h"
#include "tao/VarOut_T.h"
#include /**/ "tao/Versioned_Namespace.h"

#include "acserrC.h"

#if defined (TAO_EXPORT_MACRO)
#undef TAO_EXPORT_MACRO
#endif
#define TAO_EXPORT_MACRO 

// TAO_IDL - Generated from
// be/be_visitor_module/module_ch.cpp:49

namespace ACSErr
{
  
  // TAO_IDL - Generated from
  // be/be_visitor_constant/constant_ch.cpp:52
  
  const ACSErr::ACSErrType XBackendsErrors = 2006U;

// TAO_IDL - Generated from
// be/be_visitor_module/module_ch.cpp:78

} // module ACSErr

// TAO_IDL - Generated from
// be/be_visitor_module/module_ch.cpp:49

namespace XBackendsErrors
{
  
  // TAO_IDL - Generated from
  // be/be_visitor_constant/constant_ch.cpp:52
  
  const ACSErr::ErrorCode NoError = 0U;
  
  // TAO_IDL - Generated from
  // be/be_visitor_constant/constant_ch.cpp:52
  
  const ACSErr::ErrorCode NoSetting = 1U;
  
  // TAO_IDL - Generated from
  // be/be_visitor_exception/exception_ch.cpp:53

#if !defined (_XBACKENDSERRORS_XBACKENDSERRORSEX_CH_)
#define _XBACKENDSERRORS_XBACKENDSERRORSEX_CH_
  
  class  XBackendsErrorsEx : public ::CORBA::UserException
  {
  public:
    
    ACSErr::ErrorTrace errorTrace;
    XBackendsErrorsEx (void);
    XBackendsErrorsEx (const XBackendsErrorsEx &);
    ~XBackendsErrorsEx (void);

    XBackendsErrorsEx &operator= (const XBackendsErrorsEx &);
    
    static void _tao_any_destructor (void *);
    
    static XBackendsErrorsEx *_downcast ( ::CORBA::Exception *);
    static const XBackendsErrorsEx *_downcast ( ::CORBA::Exception const *);
    
    static ::CORBA::Exception *_alloc (void);
    
    virtual ::CORBA::Exception *_tao_duplicate (void) const;

    virtual void _raise (void) const;

    virtual void _tao_encode (TAO_OutputCDR &cdr) const;
    virtual void _tao_decode (TAO_InputCDR &cdr);
    
    
    // TAO_IDL - Generated from
    // be/be_visitor_exception/exception_ctor.cpp:66
    
    XBackendsErrorsEx (
        const ACSErr::ErrorTrace & _tao_errorTrace
      );
    
    virtual ::CORBA::TypeCode_ptr _tao_type (void) const;
  };
  
  // TAO_IDL - Generated from
  // be/be_visitor_typecode/typecode_decl.cpp:49
  
  extern  ::CORBA::TypeCode_ptr const _tc_XBackendsErrorsEx;

#endif /* end #if !defined */
  
  // TAO_IDL - Generated from
  // be/be_visitor_exception/exception_ch.cpp:53

#if !defined (_XBACKENDSERRORS_NOSETTINGEX_CH_)
#define _XBACKENDSERRORS_NOSETTINGEX_CH_
  
  class  NoSettingEx : public ::CORBA::UserException
  {
  public:
    
    ACSErr::ErrorTrace errorTrace;
    NoSettingEx (void);
    NoSettingEx (const NoSettingEx &);
    ~NoSettingEx (void);

    NoSettingEx &operator= (const NoSettingEx &);
    
    static void _tao_any_destructor (void *);
    
    static NoSettingEx *_downcast ( ::CORBA::Exception *);
    static const NoSettingEx *_downcast ( ::CORBA::Exception const *);
    
    static ::CORBA::Exception *_alloc (void);
    
    virtual ::CORBA::Exception *_tao_duplicate (void) const;

    virtual void _raise (void) const;

    virtual void _tao_encode (TAO_OutputCDR &cdr) const;
    virtual void _tao_decode (TAO_InputCDR &cdr);
    
    
    // TAO_IDL - Generated from
    // be/be_visitor_exception/exception_ctor.cpp:66
    
    NoSettingEx (
        const ACSErr::ErrorTrace & _tao_errorTrace
      );
    
    virtual ::CORBA::TypeCode_ptr _tao_type (void) const;
  };
  
  // TAO_IDL - Generated from
  // be/be_visitor_typecode/typecode_decl.cpp:49
  
  extern  ::CORBA::TypeCode_ptr const _tc_NoSettingEx;

#endif /* end #if !defined */

// TAO_IDL - Generated from
// be/be_visitor_module/module_ch.cpp:78

} // module XBackendsErrors

// TAO_IDL - Generated from
// be/be_visitor_traits.cpp:64

TAO_BEGIN_VERSIONED_NAMESPACE_DECL

// Traits specializations.
namespace TAO
{
}
TAO_END_VERSIONED_NAMESPACE_DECL



// TAO_IDL - Generated from
// be/be_visitor_exception/any_op_ch.cpp:53
TAO_BEGIN_VERSIONED_NAMESPACE_DECL



 void operator<<= (::CORBA::Any &, const XBackendsErrors::XBackendsErrorsEx &); // copying version
 void operator<<= (::CORBA::Any &, XBackendsErrors::XBackendsErrorsEx*); // noncopying version
 ::CORBA::Boolean operator>>= (const ::CORBA::Any &, XBackendsErrors::XBackendsErrorsEx *&); // deprecated
 ::CORBA::Boolean operator>>= (const ::CORBA::Any &, const XBackendsErrors::XBackendsErrorsEx *&);
TAO_END_VERSIONED_NAMESPACE_DECL



// TAO_IDL - Generated from
// be/be_visitor_exception/any_op_ch.cpp:53
TAO_BEGIN_VERSIONED_NAMESPACE_DECL



 void operator<<= (::CORBA::Any &, const XBackendsErrors::NoSettingEx &); // copying version
 void operator<<= (::CORBA::Any &, XBackendsErrors::NoSettingEx*); // noncopying version
 ::CORBA::Boolean operator>>= (const ::CORBA::Any &, XBackendsErrors::NoSettingEx *&); // deprecated
 ::CORBA::Boolean operator>>= (const ::CORBA::Any &, const XBackendsErrors::NoSettingEx *&);
TAO_END_VERSIONED_NAMESPACE_DECL



// TAO_IDL - Generated from
// be/be_visitor_exception/cdr_op_ch.cpp:52
TAO_BEGIN_VERSIONED_NAMESPACE_DECL



 ::CORBA::Boolean operator<< (TAO_OutputCDR &, const XBackendsErrors::XBackendsErrorsEx &);
 ::CORBA::Boolean operator>> (TAO_InputCDR &, XBackendsErrors::XBackendsErrorsEx &);

TAO_END_VERSIONED_NAMESPACE_DECL



// TAO_IDL - Generated from
// be/be_visitor_exception/cdr_op_ch.cpp:52
TAO_BEGIN_VERSIONED_NAMESPACE_DECL



 ::CORBA::Boolean operator<< (TAO_OutputCDR &, const XBackendsErrors::NoSettingEx &);
 ::CORBA::Boolean operator>> (TAO_InputCDR &, XBackendsErrors::NoSettingEx &);

TAO_END_VERSIONED_NAMESPACE_DECL



// TAO_IDL - Generated from
// be/be_codegen.cpp:1228
#if defined (__ACE_INLINE__)
#include "XBackendsErrorsC.inl"
#endif /* defined INLINE */

#endif /* ifndef */


