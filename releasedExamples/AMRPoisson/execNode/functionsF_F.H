#ifndef _FUNCTIONSF_F_H_
#define _FUNCTIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETPHI 
#define GUARDGETPHI 
// Prototype for Fortran procedure getphi ...
//
void FORTRAN_NAME( GETPHI ,getphi )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETPHI FORTRAN_NAME( inlineGETPHI, inlineGETPHI)
#define FORTNT_GETPHI FORTRAN_NAME( GETPHI, getphi)

inline void FORTRAN_NAME(inlineGETPHI, inlineGETPHI)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETPHI");
 FORTRAN_NAME( GETPHI ,getphi )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETPHI 

#ifndef GUARDGETMAGRESIST 
#define GUARDGETMAGRESIST 
// Prototype for Fortran procedure getmagresist ...
//
void FORTRAN_NAME( GETMAGRESIST ,getmagresist )(
      CHFp_FRA1(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) );

#define FORT_GETMAGRESIST FORTRAN_NAME( inlineGETMAGRESIST, inlineGETMAGRESIST)
#define FORTNT_GETMAGRESIST FORTRAN_NAME( GETMAGRESIST, getmagresist)

inline void FORTRAN_NAME(inlineGETMAGRESIST, inlineGETMAGRESIST)(
      CHFp_FRA1(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) )
{
 CH_TIMELEAF("FORT_GETMAGRESIST");
 FORTRAN_NAME( GETMAGRESIST ,getmagresist )(
      CHFt_FRA1(mag)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_INT(whichmag) );
}
#endif  // GUARDGETMAGRESIST 

#ifndef GUARDGETMAGPOINTRESIST 
#define GUARDGETMAGPOINTRESIST 
// Prototype for Fortran procedure getmagpointresist ...
//
void FORTRAN_NAME( GETMAGPOINTRESIST ,getmagpointresist )(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) );

#define FORT_GETMAGPOINTRESIST FORTRAN_NAME( inlineGETMAGPOINTRESIST, inlineGETMAGPOINTRESIST)
#define FORTNT_GETMAGPOINTRESIST FORTRAN_NAME( GETMAGPOINTRESIST, getmagpointresist)

inline void FORTRAN_NAME(inlineGETMAGPOINTRESIST, inlineGETMAGPOINTRESIST)(
      CHFp_REAL(mag)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_INT(icomp)
      ,CHFp_INT(whichmag) )
{
 CH_TIMELEAF("FORT_GETMAGPOINTRESIST");
 FORTRAN_NAME( GETMAGPOINTRESIST ,getmagpointresist )(
      CHFt_REAL(mag)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_INT(icomp)
      ,CHFt_INT(whichmag) );
}
#endif  // GUARDGETMAGPOINTRESIST 

#ifndef GUARDGETETARESIST 
#define GUARDGETETARESIST 
// Prototype for Fortran procedure getetaresist ...
//
void FORTRAN_NAME( GETETARESIST ,getetaresist )(
      CHFp_FRA1(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) );

#define FORT_GETETARESIST FORTRAN_NAME( inlineGETETARESIST, inlineGETETARESIST)
#define FORTNT_GETETARESIST FORTRAN_NAME( GETETARESIST, getetaresist)

inline void FORTRAN_NAME(inlineGETETARESIST, inlineGETETARESIST)(
      CHFp_FRA1(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_BOX(box)
      ,CHFp_INT(idir)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETETARESIST");
 FORTRAN_NAME( GETETARESIST ,getetaresist )(
      CHFt_FRA1(eta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_BOX(box)
      ,CHFt_INT(idir)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETETARESIST 

#ifndef GUARDGETETAPOINTRESIST 
#define GUARDGETETAPOINTRESIST 
// Prototype for Fortran procedure getetapointresist ...
//
void FORTRAN_NAME( GETETAPOINTRESIST ,getetapointresist )(
      CHFp_REAL(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_INT(idir)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) );

#define FORT_GETETAPOINTRESIST FORTRAN_NAME( inlineGETETAPOINTRESIST, inlineGETETAPOINTRESIST)
#define FORTNT_GETETAPOINTRESIST FORTRAN_NAME( GETETAPOINTRESIST, getetapointresist)

inline void FORTRAN_NAME(inlineGETETAPOINTRESIST, inlineGETETAPOINTRESIST)(
      CHFp_REAL(eta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_INT(idir)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETETAPOINTRESIST");
 FORTRAN_NAME( GETETAPOINTRESIST ,getetapointresist )(
      CHFt_REAL(eta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_INT(idir)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETETAPOINTRESIST 

#ifndef GUARDGETBETAVISCOUS 
#define GUARDGETBETAVISCOUS 
// Prototype for Fortran procedure getbetaviscous ...
//
void FORTRAN_NAME( GETBETAVISCOUS ,getbetaviscous )(
      CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(whichbeta) );

#define FORT_GETBETAVISCOUS FORTRAN_NAME( inlineGETBETAVISCOUS, inlineGETBETAVISCOUS)
#define FORTNT_GETBETAVISCOUS FORTRAN_NAME( GETBETAVISCOUS, getbetaviscous)

inline void FORTRAN_NAME(inlineGETBETAVISCOUS, inlineGETBETAVISCOUS)(
      CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_REAL(eps)
      ,CHFp_BOX(box)
      ,CHFp_INT(whichbeta) )
{
 CH_TIMELEAF("FORT_GETBETAVISCOUS");
 FORTRAN_NAME( GETBETAVISCOUS ,getbetaviscous )(
      CHFt_FRA1(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_REAL(eps)
      ,CHFt_BOX(box)
      ,CHFt_INT(whichbeta) );
}
#endif  // GUARDGETBETAVISCOUS 

#ifndef GUARDGETBETAPOINTVISCOUS 
#define GUARDGETBETAPOINTVISCOUS 
// Prototype for Fortran procedure getbetapointviscous ...
//
void FORTRAN_NAME( GETBETAPOINTVISCOUS ,getbetapointviscous )(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichbeta) );

#define FORT_GETBETAPOINTVISCOUS FORTRAN_NAME( inlineGETBETAPOINTVISCOUS, inlineGETBETAPOINTVISCOUS)
#define FORTNT_GETBETAPOINTVISCOUS FORTRAN_NAME( GETBETAPOINTVISCOUS, getbetapointviscous)

inline void FORTRAN_NAME(inlineGETBETAPOINTVISCOUS, inlineGETBETAPOINTVISCOUS)(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xval)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichbeta) )
{
 CH_TIMELEAF("FORT_GETBETAPOINTVISCOUS");
 FORTRAN_NAME( GETBETAPOINTVISCOUS ,getbetapointviscous )(
      CHFt_REAL(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xval)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichbeta) );
}
#endif  // GUARDGETBETAPOINTVISCOUS 

#ifndef GUARDGETKLBRESIST 
#define GUARDGETKLBRESIST 
// Prototype for Fortran procedure getklbresist ...
//
void FORTRAN_NAME( GETKLBRESIST ,getklbresist )(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) );

#define FORT_GETKLBRESIST FORTRAN_NAME( inlineGETKLBRESIST, inlineGETKLBRESIST)
#define FORTNT_GETKLBRESIST FORTRAN_NAME( GETKLBRESIST, getklbresist)

inline void FORTRAN_NAME(inlineGETKLBRESIST, inlineGETKLBRESIST)(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETKLBRESIST");
 FORTRAN_NAME( GETKLBRESIST ,getklbresist )(
      CHFt_FRA1(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichmag)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETKLBRESIST 

#ifndef GUARDGETKLBPOINTRESIST 
#define GUARDGETKLBPOINTRESIST 
// Prototype for Fortran procedure getklbpointresist ...
//
void FORTRAN_NAME( GETKLBPOINTRESIST ,getklbpointresist )(
      CHFp_REAL(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) );

#define FORT_GETKLBPOINTRESIST FORTRAN_NAME( inlineGETKLBPOINTRESIST, inlineGETKLBPOINTRESIST)
#define FORTNT_GETKLBPOINTRESIST FORTRAN_NAME( GETKLBPOINTRESIST, getklbpointresist)

inline void FORTRAN_NAME(inlineGETKLBPOINTRESIST, inlineGETKLBPOINTRESIST)(
      CHFp_REAL(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichmag)
      ,CHFp_INT(whicheta) )
{
 CH_TIMELEAF("FORT_GETKLBPOINTRESIST");
 FORTRAN_NAME( GETKLBPOINTRESIST ,getklbpointresist )(
      CHFt_REAL(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xvec)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichmag)
      ,CHFt_INT(whicheta) );
}
#endif  // GUARDGETKLBPOINTRESIST 

#ifndef GUARDGETKLVVISCOUS 
#define GUARDGETKLVVISCOUS 
// Prototype for Fortran procedure getklvviscous ...
//
void FORTRAN_NAME( GETKLVVISCOUS ,getklvviscous )(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_REAL(lambdafactor) );

#define FORT_GETKLVVISCOUS FORTRAN_NAME( inlineGETKLVVISCOUS, inlineGETKLVVISCOUS)
#define FORTNT_GETKLVVISCOUS FORTRAN_NAME( GETKLVVISCOUS, getklvviscous)

inline void FORTRAN_NAME(inlineGETKLVVISCOUS, inlineGETKLVVISCOUS)(
      CHFp_FRA1(klb)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_REAL(lambdafactor) )
{
 CH_TIMELEAF("FORT_GETKLVVISCOUS");
 FORTRAN_NAME( GETKLVVISCOUS ,getklvviscous )(
      CHFt_FRA1(klb)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichvel)
      ,CHFt_INT(whicheta)
      ,CHFt_INT(whichlambda)
      ,CHFt_REAL(lambdafactor) );
}
#endif  // GUARDGETKLVVISCOUS 

#ifndef GUARDGETKLVPOINTVISCOUS 
#define GUARDGETKLVPOINTVISCOUS 
// Prototype for Fortran procedure getklvpointviscous ...
//
void FORTRAN_NAME( GETKLVPOINTVISCOUS ,getklvpointviscous )(
      CHFp_REAL(klv)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_REAL(lambdafactor) );

#define FORT_GETKLVPOINTVISCOUS FORTRAN_NAME( inlineGETKLVPOINTVISCOUS, inlineGETKLVPOINTVISCOUS)
#define FORTNT_GETKLVPOINTVISCOUS FORTRAN_NAME( GETKLVPOINTVISCOUS, getklvpointviscous)

inline void FORTRAN_NAME(inlineGETKLVPOINTVISCOUS, inlineGETKLVPOINTVISCOUS)(
      CHFp_REAL(klv)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(xvec)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_INT(icomp)
      ,CHFp_REAL(eps)
      ,CHFp_INT(whichvel)
      ,CHFp_INT(whicheta)
      ,CHFp_INT(whichlambda)
      ,CHFp_REAL(lambdafactor) )
{
 CH_TIMELEAF("FORT_GETKLVPOINTVISCOUS");
 FORTRAN_NAME( GETKLVPOINTVISCOUS ,getklvpointviscous )(
      CHFt_REAL(klv)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(xvec)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_INT(icomp)
      ,CHFt_REAL(eps)
      ,CHFt_INT(whichvel)
      ,CHFt_INT(whicheta)
      ,CHFt_INT(whichlambda)
      ,CHFt_REAL(lambdafactor) );
}
#endif  // GUARDGETKLVPOINTVISCOUS 

#ifndef GUARDGETPHIPOINT 
#define GUARDGETPHIPOINT 
// Prototype for Fortran procedure getphipoint ...
//
void FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETPHIPOINT FORTRAN_NAME( inlineGETPHIPOINT, inlineGETPHIPOINT)
#define FORTNT_GETPHIPOINT FORTRAN_NAME( GETPHIPOINT, getphipoint)

inline void FORTRAN_NAME(inlineGETPHIPOINT, inlineGETPHIPOINT)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETPHIPOINT");
 FORTRAN_NAME( GETPHIPOINT ,getphipoint )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETPHIPOINT 

#ifndef GUARDGETLOFPHIZPOLY 
#define GUARDGETLOFPHIZPOLY 
// Prototype for Fortran procedure getlofphizpoly ...
//
void FORTRAN_NAME( GETLOFPHIZPOLY ,getlofphizpoly )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_GETLOFPHIZPOLY FORTRAN_NAME( inlineGETLOFPHIZPOLY, inlineGETLOFPHIZPOLY)
#define FORTNT_GETLOFPHIZPOLY FORTRAN_NAME( GETLOFPHIZPOLY, getlofphizpoly)

inline void FORTRAN_NAME(inlineGETLOFPHIZPOLY, inlineGETLOFPHIZPOLY)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_GETLOFPHIZPOLY");
 FORTRAN_NAME( GETLOFPHIZPOLY ,getlofphizpoly )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDGETLOFPHIZPOLY 

#ifndef GUARDGETPHIRZPOLY 
#define GUARDGETPHIRZPOLY 
// Prototype for Fortran procedure getphirzpoly ...
//
void FORTRAN_NAME( GETPHIRZPOLY ,getphirzpoly )(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETPHIRZPOLY FORTRAN_NAME( inlineGETPHIRZPOLY, inlineGETPHIRZPOLY)
#define FORTNT_GETPHIRZPOLY FORTRAN_NAME( GETPHIRZPOLY, getphirzpoly)

inline void FORTRAN_NAME(inlineGETPHIRZPOLY, inlineGETPHIRZPOLY)(
      CHFp_REAL(phi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETPHIRZPOLY");
 FORTRAN_NAME( GETPHIRZPOLY ,getphirzpoly )(
      CHFt_REAL(phi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETPHIRZPOLY 

#ifndef GUARDGETGRADPHIRZPOLY 
#define GUARDGETGRADPHIRZPOLY 
// Prototype for Fortran procedure getgradphirzpoly ...
//
void FORTRAN_NAME( GETGRADPHIRZPOLY ,getgradphirzpoly )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADPHIRZPOLY FORTRAN_NAME( inlineGETGRADPHIRZPOLY, inlineGETGRADPHIRZPOLY)
#define FORTNT_GETGRADPHIRZPOLY FORTRAN_NAME( GETGRADPHIRZPOLY, getgradphirzpoly)

inline void FORTRAN_NAME(inlineGETGRADPHIRZPOLY, inlineGETGRADPHIRZPOLY)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADPHIRZPOLY");
 FORTRAN_NAME( GETGRADPHIRZPOLY ,getgradphirzpoly )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADPHIRZPOLY 

#ifndef GUARDGETGRADPHIPOINT 
#define GUARDGETGRADPHIPOINT 
// Prototype for Fortran procedure getgradphipoint ...
//
void FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADPHIPOINT FORTRAN_NAME( inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)
#define FORTNT_GETGRADPHIPOINT FORTRAN_NAME( GETGRADPHIPOINT, getgradphipoint)

inline void FORTRAN_NAME(inlineGETGRADPHIPOINT, inlineGETGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADPHIPOINT");
 FORTRAN_NAME( GETGRADPHIPOINT ,getgradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADPHIPOINT 

#ifndef GUARDGETLOFPHI 
#define GUARDGETLOFPHI 
// Prototype for Fortran procedure getlofphi ...
//
void FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) );

#define FORT_GETLOFPHI FORTRAN_NAME( inlineGETLOFPHI, inlineGETLOFPHI)
#define FORTNT_GETLOFPHI FORTRAN_NAME( GETLOFPHI, getlofphi)

inline void FORTRAN_NAME(inlineGETLOFPHI, inlineGETLOFPHI)(
      CHFp_FRA1(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETLOFPHI");
 FORTRAN_NAME( GETLOFPHI ,getlofphi )(
      CHFt_FRA1(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETLOFPHI 

#ifndef GUARDGETLOFPHIPOINT 
#define GUARDGETLOFPHIPOINT 
// Prototype for Fortran procedure getlofphipoint ...
//
void FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_GETLOFPHIPOINT FORTRAN_NAME( inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)
#define FORTNT_GETLOFPHIPOINT FORTRAN_NAME( GETLOFPHIPOINT, getlofphipoint)

inline void FORTRAN_NAME(inlineGETLOFPHIPOINT, inlineGETLOFPHIPOINT)(
      CHFp_REAL(lofphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_GETLOFPHIPOINT");
 FORTRAN_NAME( GETLOFPHIPOINT ,getlofphipoint )(
      CHFt_REAL(lofphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDGETLOFPHIPOINT 

#ifndef GUARDGETDBGPHI 
#define GUARDGETDBGPHI 
// Prototype for Fortran procedure getdbgphi ...
//
void FORTRAN_NAME( GETDBGPHI ,getdbgphi )(
      CHFp_FRA1(dbgphi)
      ,CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_BOX(box) );

#define FORT_GETDBGPHI FORTRAN_NAME( inlineGETDBGPHI, inlineGETDBGPHI)
#define FORTNT_GETDBGPHI FORTRAN_NAME( GETDBGPHI, getdbgphi)

inline void FORTRAN_NAME(inlineGETDBGPHI, inlineGETDBGPHI)(
      CHFp_FRA1(dbgphi)
      ,CHFp_FRA1(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETDBGPHI");
 FORTRAN_NAME( GETDBGPHI ,getdbgphi )(
      CHFt_FRA1(dbgphi)
      ,CHFt_FRA1(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETDBGPHI 

#ifndef GUARDGETDBGPHIPOINT 
#define GUARDGETDBGPHIPOINT 
// Prototype for Fortran procedure getdbgphipoint ...
//
void FORTRAN_NAME( GETDBGPHIPOINT ,getdbgphipoint )(
      CHFp_REAL(dbgphi)
      ,CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha) );

#define FORT_GETDBGPHIPOINT FORTRAN_NAME( inlineGETDBGPHIPOINT, inlineGETDBGPHIPOINT)
#define FORTNT_GETDBGPHIPOINT FORTRAN_NAME( GETDBGPHIPOINT, getdbgphipoint)

inline void FORTRAN_NAME(inlineGETDBGPHIPOINT, inlineGETDBGPHIPOINT)(
      CHFp_REAL(dbgphi)
      ,CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(alpha) )
{
 CH_TIMELEAF("FORT_GETDBGPHIPOINT");
 FORTRAN_NAME( GETDBGPHIPOINT ,getdbgphipoint )(
      CHFt_REAL(dbgphi)
      ,CHFt_REAL(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(alpha) );
}
#endif  // GUARDGETDBGPHIPOINT 

#ifndef GUARDGETBETAPOINT 
#define GUARDGETBETAPOINT 
// Prototype for Fortran procedure getbetapoint ...
//
void FORTRAN_NAME( GETBETAPOINT ,getbetapoint )(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETBETAPOINT FORTRAN_NAME( inlineGETBETAPOINT, inlineGETBETAPOINT)
#define FORTNT_GETBETAPOINT FORTRAN_NAME( GETBETAPOINT, getbetapoint)

inline void FORTRAN_NAME(inlineGETBETAPOINT, inlineGETBETAPOINT)(
      CHFp_REAL(beta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETBETAPOINT");
 FORTRAN_NAME( GETBETAPOINT ,getbetapoint )(
      CHFt_REAL(beta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETBETAPOINT 

#ifndef GUARDGETGRADBETAPOINT 
#define GUARDGETGRADBETAPOINT 
// Prototype for Fortran procedure getgradbetapoint ...
//
void FORTRAN_NAME( GETGRADBETAPOINT ,getgradbetapoint )(
      CHFp_REALVECT(gradbeta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETGRADBETAPOINT FORTRAN_NAME( inlineGETGRADBETAPOINT, inlineGETGRADBETAPOINT)
#define FORTNT_GETGRADBETAPOINT FORTRAN_NAME( GETGRADBETAPOINT, getgradbetapoint)

inline void FORTRAN_NAME(inlineGETGRADBETAPOINT, inlineGETGRADBETAPOINT)(
      CHFp_REALVECT(gradbeta)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETGRADBETAPOINT");
 FORTRAN_NAME( GETGRADBETAPOINT ,getgradbetapoint )(
      CHFt_REALVECT(gradbeta)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETGRADBETAPOINT 

#ifndef GUARDGETBETAGRADPHIPOINT 
#define GUARDGETBETAGRADPHIPOINT 
// Prototype for Fortran procedure getbetagradphipoint ...
//
void FORTRAN_NAME( GETBETAGRADPHIPOINT ,getbetagradphipoint )(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_GETBETAGRADPHIPOINT FORTRAN_NAME( inlineGETBETAGRADPHIPOINT, inlineGETBETAGRADPHIPOINT)
#define FORTNT_GETBETAGRADPHIPOINT FORTRAN_NAME( GETBETAGRADPHIPOINT, getbetagradphipoint)

inline void FORTRAN_NAME(inlineGETBETAGRADPHIPOINT, inlineGETBETAGRADPHIPOINT)(
      CHFp_REALVECT(gradphi)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_GETBETAGRADPHIPOINT");
 FORTRAN_NAME( GETBETAGRADPHIPOINT ,getbetagradphipoint )(
      CHFt_REALVECT(gradphi)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDGETBETAGRADPHIPOINT 

#ifndef GUARDGETSRC 
#define GUARDGETSRC 
// Prototype for Fortran procedure getsrc ...
//
void FORTRAN_NAME( GETSRC ,getsrc )(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) );

#define FORT_GETSRC FORTRAN_NAME( inlineGETSRC, inlineGETSRC)
#define FORTNT_GETSRC FORTRAN_NAME( GETSRC, getsrc)

inline void FORTRAN_NAME(inlineGETSRC, inlineGETSRC)(
      CHFp_FRA1(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(diffconst)
      ,CHFp_CONST_REALVECT(problo)
      ,CHFp_CONST_REALVECT(probhi)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETSRC");
 FORTRAN_NAME( GETSRC ,getsrc )(
      CHFt_FRA1(src)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(diffconst)
      ,CHFt_CONST_REALVECT(problo)
      ,CHFt_CONST_REALVECT(probhi)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETSRC 

#ifndef GUARDGETSRCPOINT 
#define GUARDGETSRCPOINT 
// Prototype for Fortran procedure getsrcpoint ...
//
void FORTRAN_NAME( GETSRCPOINT ,getsrcpoint )(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(diffconst) );

#define FORT_GETSRCPOINT FORTRAN_NAME( inlineGETSRCPOINT, inlineGETSRCPOINT)
#define FORTNT_GETSRCPOINT FORTRAN_NAME( GETSRCPOINT, getsrcpoint)

inline void FORTRAN_NAME(inlineGETSRCPOINT, inlineGETSRCPOINT)(
      CHFp_REAL(src)
      ,CHFp_CONST_REALVECT(freq)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(diffconst) )
{
 CH_TIMELEAF("FORT_GETSRCPOINT");
 FORTRAN_NAME( GETSRCPOINT ,getsrcpoint )(
      CHFt_REAL(src)
      ,CHFt_CONST_REALVECT(freq)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(diffconst) );
}
#endif  // GUARDGETSRCPOINT 

}

#endif
