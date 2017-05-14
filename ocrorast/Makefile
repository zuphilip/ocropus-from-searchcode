CXXFLAGS=-fPIC -g -O4

libocrorast.a: \
	wrappers.o line-info.o linesegs.o log-reg-data.o ocr-char-stats.o \
	ocr-classify-zones.o ocr-color-encode-layout.o \
        ocr-ctextline-rast-extended.o \
	ocr-ctextline-rast.o ocr-deskew-rast.o \
        ocr-detect-columns.o ocr-detect-paragraphs.o \
	ocr-doc-clean-concomp.o ocr-doc-clean.o \
        ocr-extract-gutters.o ocr-extract-rulings.o \
	ocr-layout-1cp.o ocr-layout-rast.o \
        ocr-layout-rast1.o ocr-layout-smear.o \
	ocr-noisefilter.o ocr-pageframe-rast.o \
        ocr-pageseg-wcuts.o ocr-pageseg-xycut.o \
	ocr-reading-order.o ocr-segmentations.o \
        ocr-text-image-seg.o ocr-utils.o \
	ocr-visualize-layout-rast.o ocr-whitespace-cover.o \
        ocr-word-segmentation.o pagesegs.o
	ar cr $@ $^

clean:
	rm -rf build *.o *.so *.a *_wrap.cpp ocrorast.py
