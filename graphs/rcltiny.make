
rcl:
	rcl.sh setup tiny -n rcltiny.mci -t rcltiny.tab
	rcl.sh tree tiny rcltiny.cls?
	RCL_RES_DOT_TREE=16 rcl.sh res tiny -r "1 2 3 4 5"

tree: tiny/rcl.joindot
	RCL_DOT_TREE_ANNOT=theannot.txt RCL_DOT_TREE_HEADER=theheader.txt dot-rcl-mergetree.pl tiny/rcl.joindot > tiny/j.dot &&  dot -Tpdf -Gsize=5,5\! < tiny/j.dot > tiny/j.pdf
	gs -o tiny-tree.pdf -sDEVICE=pdfwrite -dColorConversionStrategy=/sRGB -dProcessColorModel=/DeviceRGB tiny/j.pdf

grid:
	mkdir -p tiny
	dot -Tpdf -Gsize=5,5\! < rcltiny1.netdot > tiny/c1.pdf
	dot -Tpdf -Gsize=5,5\! < rcltiny2.netdot > tiny/c2.pdf
	dot -Tpdf -Gsize=5,5\! < rcltiny3.netdot > tiny/c3.pdf
	gs -o tiny-cls1.pdf -sDEVICE=pdfwrite -dColorConversionStrategy=/sRGB -dProcessColorModel=/DeviceRGB tiny/c1.pdf
	gs -o tiny-cls2.pdf -sDEVICE=pdfwrite -dColorConversionStrategy=/sRGB -dProcessColorModel=/DeviceRGB tiny/c2.pdf
	gs -o tiny-cls3.pdf -sDEVICE=pdfwrite -dColorConversionStrategy=/sRGB -dProcessColorModel=/DeviceRGB tiny/c3.pdf

