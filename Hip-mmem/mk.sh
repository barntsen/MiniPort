#!/bin/sh
opt=

hipcc -O3  $opt -c ac2dmod.cpp
hipcc -O3  $opt -c ac2d.cpp
hipcc -O3  $opt -c diff.cpp
hipcc -O3  $opt -c model.cpp
hipcc -O3  $opt -c rec.cpp
hipcc -O3  $opt -c src.cpp
hipcc -O3  $opt -c util.cpp
hipcc -O3  $opt -c gpu.cpp

ar cr libac2d.a ac2d.o diff.o model.o rec.o src.o util.o gpu.o

hipcc  -o ac2dmod ac2dmod.o ac2d.o model.o rec.o src.o util.o diff.o gpu.o 
