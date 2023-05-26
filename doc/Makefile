NAME = manual
TYPE = pdf
OUTPUT = $(NAME).$(TYPE)

all: $(OUTPUT)

%.pdf: %.texi
	texi2pdf $<

clean:
	rm -f *.log *.aux *.cp *.fn *.ky *.tp *.pg *.vr *.toc

veryclean: clean
	rm -f $(OUTPUT) 

