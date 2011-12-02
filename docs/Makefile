all: reference.html 

%.png: %.svg
	rsvg $< $@

%.html: %.rst
	rst2html --stylesheet=lsr.css --link-stylesheet $< $@
