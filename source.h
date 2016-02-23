######################################################
#user should specify their own source code here      #
#user should not changed the macros(those in upcase) #
#the major Makefile in bin/                          #
#the two supporting Makefile are in fobj/ and cobj/  #
#user need to put them in the specified dir          #
#if user's the dir structures are the same as the    #
# example, they don't need to change the Makefiles   #
# only need to change source.h                       #
######################################################

CSRC1 = delete.C depchk.C \
        hadpt.C  \
        updatenei.C update_info_interp.C \
        unrefine.C
CSRC2 = element2.C hashtab2.C node.C
CSRC3 = assemble.C \
	assemble_bubble.C 
CSRC4 = compare_key.C datread.C delete_tab.C \
	hpfem.C new_datatype.C hilbert.C refine2.C 
CSRC5 = #bb_bb.C gldof.C nn_ss.C ordering.C rmdof.C vv_ee.C
CSRC6 = GisApi.C GisBinFile.C GisAscFile.C GisRasterHdr.C
CSRC7 = fill.C update_element_info.C 
CSRC8 = tecplot.C \
        BSFC_combine_elements.C repartition_BSFC.C \
        BSFC_create_bins.C BSFC_refine_partition.C \
	BSFC_update_and_send_elements.C \
	BSFC_update_element_proc.C
CSRC9 = 
CSRC10 = step.C calc_volume.C \
         move_data.C slopelimit.C \
         element_weight.C
#cSRC1  = GrassApi.c

CSRC  = $(CSRC1) $(CSRC2) $(CSRC3) $(CSRC4) $(CSRC5) \
	$(CSRC6) $(CSRC7) $(CSRC8) $(CSRC9) $(CSRC10) \
	$(CSRC11) 
cSRC  = $(cSRC1)

FSRC1 = bcmax.f bcoeff.f shapeb.f
FSRC2 = elmcon.f mcoeff.f  \
        masselem.f ela_coef.f
FSRC3 = maxval.f minval.f
FSRC4 = bointg.f cnst.f setz.f setzi.f
FSRC5 = aslmv2.f rhsub.f sccallp.f schur5.f \
	transf.f tri.f triev.f eigsrt.f
FSRC6 = eval2.f 
FSRC7 = dshap2dg.f gshape.f shape2dg.f shape1.f \
        shape2.f dshap1.f dshap2.f setshape.f
FSRC8 =
FSRC9 = gauss2.f gausse.f decomp.f 
FSRC10 = exsol.f  
FSRC11 = getcoef.f sgn.f eigen.f predict.f correct.f \
	 elbshal.f getquad.f getkactxy.f

FSRC  = $(FSRC1) $(FSRC2) $(FSRC3) $(FSRC4) $(FSRC5) \
	$(FSRC6) $(FSRC7) $(FSRC8) $(FSRC9) $(FSRC10) \
        $(FSRC11)

CSRCDIR = ../csrc/
FSRCDIR = ../fsrc/

CUSERDIR = adapt/ datstr/ getsol/ main/ ordering/ \
	   postproc/ repartition/ tecplot/ block/ \
           geoflow/ 
FUSERDIR = pack_bc/ pack_elemn/ pack_err/ pack_initial/ \
	   pack_locit/ pack_post/ pack_shap/ pack_stress/ \
	   pack_util/ pack_custom/ pack_geoflow/

COBJDIR = ../cobj/
FOBJDIR = ../fobj/
cOBJDIR = ../ccobj/

##### end user defination #####
