all: Ellip3dReaderGUI.bqrc

Ellip3dReaderGUI.bqrc: gui.qrc Ellip3dReader.ui fileTypes.xml
	rcc -binary -o Ellip3dReaderGUI.bqrc gui.qrc

clean:
	rm -f Ellip3dReaderGUI.bqrc
