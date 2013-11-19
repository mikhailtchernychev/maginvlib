# Microsoft Developer Studio Project File - Name="maginvlib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=maginvlib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "maginvlib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "maginvlib.mak" CFG="maginvlib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "maginvlib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "maginvlib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "maginvlib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /Zp1 /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "maginvlib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /Zp1 /MTd /W3 /Gm /GR /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "maginvlib - Win32 Release"
# Name "maginvlib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\dcov.cpp
# End Source File
# Begin Source File

SOURCE=.\df551_1.c
# End Source File
# Begin Source File

SOURCE=.\dipole1d.cpp
# End Source File
# Begin Source File

SOURCE=.\dnls1.cpp
# End Source File
# Begin Source File

SOURCE=.\dnls1estimator.cpp
# End Source File
# Begin Source File

SOURCE=.\dnls1m_dr.cpp
# End Source File
# Begin Source File

SOURCE=.\dz29rls.cpp
# End Source File
# Begin Source File

SOURCE=.\dz29svd.cpp
# End Source File
# Begin Source File

SOURCE=.\gradprofile.cpp
# End Source File
# Begin Source File

SOURCE=.\hooke.c
# End Source File
# Begin Source File

SOURCE=.\invers_m.c
# End Source File
# Begin Source File

SOURCE=.\l1estimator.cpp
# End Source File
# Begin Source File

SOURCE=.\l2estimator.cpp
# End Source File
# Begin Source File

SOURCE=.\line.cpp
# End Source File
# Begin Source File

SOURCE=.\magarray.cpp
# End Source File
# Begin Source File

SOURCE=.\magarrayt.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\magbackground.cpp
# End Source File
# Begin Source File

SOURCE=.\magbackground2.cpp
# End Source File
# Begin Source File

SOURCE=.\magchannel.cpp
# End Source File
# Begin Source File

SOURCE=.\magcollection.cpp
# End Source File
# Begin Source File

SOURCE=.\magdata.cpp
# End Source File
# Begin Source File

SOURCE=.\magdata1d.cpp
# End Source File
# Begin Source File

SOURCE=.\magdatacollection.cpp
# End Source File
# Begin Source File

SOURCE=.\magdipole.cpp
# End Source File
# Begin Source File

SOURCE=.\magdipolethread.cpp
# End Source File
# Begin Source File

SOURCE=.\magdipolethread.h
# End Source File
# Begin Source File

SOURCE=.\magfunc.cpp
# End Source File
# Begin Source File

SOURCE=.\maggriddata.cpp
# End Source File
# Begin Source File

SOURCE=.\maginv.cpp
# End Source File
# Begin Source File

SOURCE=.\magobject.cpp
# End Source File
# Begin Source File

SOURCE=.\magpipe.cpp
# End Source File
# Begin Source File

SOURCE=.\magprofiles.cpp
# End Source File
# Begin Source File

SOURCE=.\magsegments.cpp
# End Source File
# Begin Source File

SOURCE=.\magsnake.cpp
# End Source File
# Begin Source File

SOURCE=.\magtest.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\modellines.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\monopole.cpp
# End Source File
# Begin Source File

SOURCE=.\splik.cpp
# End Source File
# Begin Source File

SOURCE=.\testprofile.cpp
# End Source File
# Begin Source File

SOURCE=.\tn_dr.cpp
# End Source File
# Begin Source File

SOURCE=.\tn_estimator.cpp
# End Source File
# Begin Source File

SOURCE=.\tn_m.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\arrayinv.h
# End Source File
# Begin Source File

SOURCE=.\dipole1d.h
# End Source File
# Begin Source File

SOURCE=.\dnls1estimator.h
# End Source File
# Begin Source File

SOURCE=.\dnls1m_dr.h
# End Source File
# Begin Source File

SOURCE=.\fcnbase.h
# End Source File
# Begin Source File

SOURCE=.\fieldtype.h
# End Source File
# Begin Source File

SOURCE=.\gradprofile.h
# End Source File
# Begin Source File

SOURCE=.\grd7.h
# End Source File
# Begin Source File

SOURCE=.\inversion.h
# End Source File
# Begin Source File

SOURCE=.\l1estimator.h
# End Source File
# Begin Source File

SOURCE=.\l2estimator.h
# End Source File
# Begin Source File

SOURCE=.\magarray.h
# End Source File
# Begin Source File

SOURCE=.\magarrayt.h
# End Source File
# Begin Source File

SOURCE=.\magbackground.h
# End Source File
# Begin Source File

SOURCE=.\magbackground2.h
# End Source File
# Begin Source File

SOURCE=.\magchannel.h
# End Source File
# Begin Source File

SOURCE=.\magcollection.h
# End Source File
# Begin Source File

SOURCE=.\magdata.h
# End Source File
# Begin Source File

SOURCE=.\magdata1d.h
# End Source File
# Begin Source File

SOURCE=.\magdata1dt.h
# End Source File
# Begin Source File

SOURCE=.\magdatacollection.h
# End Source File
# Begin Source File

SOURCE=.\magdipole.h
# End Source File
# Begin Source File

SOURCE=.\magfunc.h
# End Source File
# Begin Source File

SOURCE=.\maggriddata.h
# End Source File
# Begin Source File

SOURCE=.\maginv.h
# End Source File
# Begin Source File

SOURCE=.\magobject.h
# End Source File
# Begin Source File

SOURCE=.\magpipe.h
# End Source File
# Begin Source File

SOURCE=.\magprofiles.h
# End Source File
# Begin Source File

SOURCE=.\magsegments.h
# End Source File
# Begin Source File

SOURCE=.\magsnake.h
# End Source File
# Begin Source File

SOURCE=.\monopole.h
# End Source File
# Begin Source File

SOURCE=.\profwrapper.h
# End Source File
# Begin Source File

SOURCE=.\testprofile.h
# End Source File
# Begin Source File

SOURCE=.\tn_dr.h
# End Source File
# Begin Source File

SOURCE=.\tn_estimator.h
# End Source File
# End Group
# End Target
# End Project
