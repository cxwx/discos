#*******************************************************************************
# PPPPPPPP
#
# "@(#) $Id: Makefile,v 1.2 2007-06-05 06:26:22 c.migoni Exp $"
#
# Makefile of ........
#
# who       when      what
# --------  --------  ----------------------------------------------
# almamgr  28/07/05  created
#

#*******************************************************************************
# This Makefile follows VLT Standards (see Makefile(5) for more).
#*******************************************************************************
# REMARKS
#    None
#------------------------------------------------------------------------

#
# user definable C-compilation flags
USER_CFLAGS = -c -pedantic -Wall

#
# additional include and library search paths
#USER_INC = 
USER_LIB = -lm 

#
# MODULE CODE DESCRIPTION:
# ------------------------
# As a general rule:  public file are "cleaned" and "installed"  
#                     local (_L) are not "installed".

#
# C programs (public and local)
# -----------------------------
EXECUTABLES     =
EXECUTABLES_L   = 

#
# <brief description of xxxxx program>
xxxxx_OBJECTS   =	
xxxxx_LDFLAGS   =
xxxxx_LIBS      =

#
# special compilation flags for single c sources
#yyyyy_CFLAGS   = 

#
# Includes (.h) files (public only)
# ---------------------------------
INCLUDES        = slalib.h slamac.h

#
# Libraries (public and local)
# ----------------------------
LIBRARIES       = SlaLibrary
LIBRARIES_L     =
SlaLibrary_OBJECTS = addet afin airmas altaz amp \
	ampqk aop aoppa aoppat aopqk atmdsp \
	av2m bear c2i c2ipa c2ipad c2ipas c2iqk caf2r c2iqkz caldj calyd cc2s \
	cc62s cd2tf cldj clyd combn cr2af cr2tf \
	cs2c cs2c6 ctf2d ctf2r daf2r dafin \
	dat dav2m dbear dbjin dc62s dcc2s \
	dcmpf dcs2c dd2tf de2h deuler dfltin \
	dh2e dimxv djcal djcl dm2av dmat \
	dmoon dmxm dmxv dpav dr2af dr2tf \
	drange dranrm ds2c6 ds2tp dsep dsepv dt \
	dtf2d dtf2r dtp2s dtp2v dtps2c dtpv2c \
	dtt dv2tp dvdv dvn dvxv e2h earth \
	ecleq ecmat ecor eg50 el2ue eo eors epb \
	epb2d epco epj epj2d epv eqecl eqeqx \
	eqgal era etrms euler evp fitxy fk425 \
	fk45z fk524 fk52h fk54z fk5hz flotin fw2m fw2xy \
	g2ixys galeq galsup ge50 geoc gmst gmsta gst \
	h2e h2fk5 hfk5z imxv intin invf \
	i2c i2cqk i2o i2opa i2opad i2opat i2oqk int2in \
	kbj m2av map mappa mapqk mapqkz \
	moon mxm mxv nut nutc oap oapqk \
	nu00a nu nutc80\
	o2i o2iqk obs pa pav pcd pda2h pdq2h permut \
	pertel pertue pfw planel planet plante plantu pm pncio pneqx\
	polmo pomom prebn prec precl preces prenut \
	pvobs pv2el pv2ue pxy range ranorm \
	rcc rdplan refco refcoq refro refv \
	refz rverot rvgalc rvlg rvlsrd rvlsrk \
	s s2tp sep sepv smat sp subet supgal \
	svd svdcov svdsol tp2s tp2v tps2c \
	tpv2c ue2el ue2pv unpcd v2tp vdv \
	vn vxv xy2xy zd 
SlaLibrary_LIBS =

#
# <brief description of lllll library>
lllll_OBJECTS   =

#
# Scripts (public and local)
# ----------------------------
SCRIPTS         =
SCRIPTS_L       =

#
# TCL scripts (public and local)
# ------------------------------
TCL_SCRIPTS     =
TCL_SCRIPTS_L   =

#
# Python stuff (public and local)
# ----------------------------
PY_SCRIPTS         =
PY_SCRIPTS_L       =

PY_MODULES         =
PY_MODULES_L       =

PY_PACKAGES        =
PY_PACKAGES_L      =
pppppp_MODULES	   =

#
# <brief description of tttttt tcl-script>
tttttt_OBJECTS  =
tttttt_TCLSH    = 
tttttt_LIBS     = 

#
# TCL libraries (public and local)
# ------------------------------
TCL_LIBRARIES   =
TCL_LIBRARIES_L =

#
# <brief description of tttlll library>
tttlll_OBJECTS  = 

#
# Configuration Database Files
# ----------------------------
CDB_SCHEMAS = 

# 
# IDL Files and flags
# 
IDL_FILES =
IDL_TAO_FLAGS =
USER_IDL =
#
# Jarfiles and their directories
#
JARFILES= 
jjj_DIRS=
jjj_EXTRAS= 
#
# java sources in Jarfile on/off
DEBUG= 
#
# ACS XmlIdl generation on/off
#
XML_IDL= 
#
# Java Component Helper Classes generation on/off
#
COMPONENT_HELPERS=
#
# Java Entity Classes generation on/off
#
XSDBIND=
#
# Schema Config files for the above
#
XSDBIND_INCLUDE=
# man pages to be done
# --------------------
MANSECTIONS =
MAN1 =
MAN3 =
MAN5 =
MAN7 =
MAN8 =

#
# local man pages
# ---------------
MANl =

#
# ASCII file to be converted into Framemaker-MIF
# --------------------
ASCII_TO_MIF = 

#
# other files to be installed
#----------------------------
INSTALL_FILES =

#
# list of all possible C-sources (used to create automatic dependencies)
# ------------------------------
CSOURCENAMES = \
	$(foreach exe, $(EXECUTABLES) $(EXECUTABLES_L), $($(exe)_OBJECTS)) \
	$(foreach rtos, $(RTAI_MODULES) , $($(rtos)_OBJECTS)) \
	$(foreach lib, $(LIBRARIES) $(LIBRARIES_L), $($(lib)_OBJECTS))

#
#>>>>> END OF standard rules

#
# INCLUDE STANDARDS
# -----------------

MAKEDIRTMP := $(shell searchFile include/acsMakefile)
ifneq ($(MAKEDIRTMP),\#error\#)
   MAKEDIR := $(MAKEDIRTMP)/include
   include $(MAKEDIR)/acsMakefile
endif

#
# TARGETS
# -------
all:	do_all
	@echo " . . . 'all' done" 

clean : clean_all 
	@echo " . . . clean done"

clean_dist : clean_all clean_dist_all 
	@echo " . . . clean_dist done"

man   : do_man 
	@echo " . . . man page(s) done"

install : install_all
	@echo " . . . installation done"


#___oOo___
