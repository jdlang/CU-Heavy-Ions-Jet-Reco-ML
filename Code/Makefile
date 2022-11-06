# Add here if multiple files need to be compiled
all:
	jordanevents
	joeyevents

# Compiles Event_Generator_PYTHIA.C regardless of changes
jordanevents: Event_Generator_PYTHIA.C
	g++ Event_Generator_PYTHIA.C /home/qq/Documents/pythia/pythia8307/lib/libpythia8.a -o Event_Generator_PYTHIA -I/home/qq/fastjet-install/include/ -I/home/qq/Documents/pythia/pythia8307/include/  -L/home/qq/fastjet-install/lib/ -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,/home/qq/Documents/pythia/pythia8307/lib/ -ldl `/home/qq/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -ldl `root-config --libs --cflags`

# NOTE: The above paths need to be updated for local pythia8 installation and fastjet installation

joeyevents: joey_eventgen.C
	g++ joey_eventgen.C /home/qq/Documents/pythia/pythia8307/lib/libpythia8.a -o Event_Generator_PYTHIA -I/home/qq/fastjet-install/include/ -I/home/qq/Documents/pythia/pythia8307/include/  -L/home/qq/fastjet-install/lib/ -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath,/home/qq/Documents/pythia/pythia8307/lib/ -ldl `/home/qq/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -ldl `root-config --libs --cflags`

