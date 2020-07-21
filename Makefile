

OBJECTS = qpalloc.o qpdifmat0.o   qpsmat0.o    qpstart4.o   spbdphj.o \
banddpasss.o   qpdifmatl.o   qpsmatc.o    qpstart4a.o   spbdphy.o \
butterworth.o  qpdifmats.o   qpsprop.o    qpstart4g.o   spbh12.o \
caxcb.o        qpfftinv.o    qpsprop0.o   qpstart6.o    spbjh.o \
cdsvd500.o     qpgetinp.o    qpsprop1.o   qpstart6a.o   spbphj.o \
cmemcpy.o      qpgrnspec.o   qpspropg.o   qpstart6g.o   spbphy.o \
disazi.o        qpspropg0.o  qpsublayer.o  spbpsj.o \
four1w.o       qppsvkern.o   qpspropg1.o  qptmat.o      spbpsy.o \
fsimpson.o     qppsvkerng.o  qpstart0.o   qptprop.o     swavelet.o \
legendre.o     qpqmodel.o    qpstart0a.o  qpwvint.o     taper.o \
moments.o      qpshkern.o    qpstart0g.o  ruku.o        transfs2t.o \
      qpsmat.o      qpstart2t.o  skipdoc.o


install: qssp 
	mv qssp ../bin

qssp: $(OBJECTS) qpmain.o
	gfortran  $(OBJECTS) qpmain.o -o qssp

qpmain.o : qpmain.f90
	gfortran -c $<

%.o : %.f
	gfortran  -c  $< 

clean:
	rm -f *.o *.mod