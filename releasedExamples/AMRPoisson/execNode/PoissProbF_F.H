#ifndef _POISSPROBF_F_H_
#define _POISSPROBF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETRHSPOIS 
#define GUARDGETRHSPOIS 
// Prototype for Fortran procedure GETRHSPOIS ...
//
void FORTRAN_NAME( GETRHSPOIS ,getrhspois )(
      CHFp_FRA(rhs)
      ,CHFp_BOX(box)
      ,CHFp_BOX(domain)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(rhono)
      ,CHFp_CONST_REAL(rno)
      ,CHFp_CONST_INT(iprob) );

#define FORT_GETRHSPOIS FORTRAN_NAME( inlineGETRHSPOIS, inlineGETRHSPOIS)
#define FORTNT_GETRHSPOIS FORTRAN_NAME( GETRHSPOIS, getrhspois)

inline void FORTRAN_NAME(inlineGETRHSPOIS, inlineGETRHSPOIS)(
      CHFp_FRA(rhs)
      ,CHFp_BOX(box)
      ,CHFp_BOX(domain)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(rhono)
      ,CHFp_CONST_REAL(rno)
      ,CHFp_CONST_INT(iprob) )
{
 CH_TIMELEAF("FORT_GETRHSPOIS");
 FORTRAN_NAME( GETRHSPOIS ,getrhspois )(
      CHFt_FRA(rhs)
      ,CHFt_BOX(box)
      ,CHFt_BOX(domain)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(rhono)
      ,CHFt_CONST_REAL(rno)
      ,CHFt_CONST_INT(iprob) );
}
#endif  // GUARDGETRHSPOIS 

#ifndef GUARDGETRHSEXACT 
#define GUARDGETRHSEXACT 
// Prototype for Fortran procedure getrhsexact ...
//
void FORTRAN_NAME( GETRHSEXACT ,getrhsexact )(
      CHFp_REAL(rhono)
      ,CHFp_CONST_REAL(rno)
      ,CHFp_CONST_REAL(rad)
      ,CHFp_CONST_INT(iprob)
      ,CHFp_REAL(rhsexact) );

#define FORT_GETRHSEXACT FORTRAN_NAME( inlineGETRHSEXACT, inlineGETRHSEXACT)
#define FORTNT_GETRHSEXACT FORTRAN_NAME( GETRHSEXACT, getrhsexact)

inline void FORTRAN_NAME(inlineGETRHSEXACT, inlineGETRHSEXACT)(
      CHFp_REAL(rhono)
      ,CHFp_CONST_REAL(rno)
      ,CHFp_CONST_REAL(rad)
      ,CHFp_CONST_INT(iprob)
      ,CHFp_REAL(rhsexact) )
{
 CH_TIMELEAF("FORT_GETRHSEXACT");
 FORTRAN_NAME( GETRHSEXACT ,getrhsexact )(
      CHFt_REAL(rhono)
      ,CHFt_CONST_REAL(rno)
      ,CHFt_CONST_REAL(rad)
      ,CHFt_CONST_INT(iprob)
      ,CHFt_REAL(rhsexact) );
}
#endif  // GUARDGETRHSEXACT 

}

#endif
