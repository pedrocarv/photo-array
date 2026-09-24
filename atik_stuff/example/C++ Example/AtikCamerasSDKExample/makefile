SRC := AtikCamerasSDKApp.cpp ConsoleHelper.cpp main.cpp

INCLUDE_DIR = -I./

CFLAGS = $(INCLUDE_DIR) -Wall -Wextra -std=c++14

LDFLAGS = -latikcameras -lusb-1.0 -ludev
#In case user wants to keep the libatikcameras.so locally, then copy the corresponding libatikcameras.so to this folder, uncomment the LDFLAGS below and comment the LDFLAGS from above
#LDFLAGS = -Wl,-rpath=./ -L./ -latikcameras -lusb-1.0 -ludev

OBJ_DIR := obj/
OBJ := $(patsubst %.cpp,$(OBJ_DIR)%.o, $(SRC))

BUILDDIR := bin/

$(BUILDDIR)atiksdkexample.out: $(OBJ)
	@mkdir -p $(BUILDDIR)
	$(CXX) -o $@ $^ $(LDFLAGS) 

$(OBJ_DIR)%.o: %.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) -o $@ -c $< $(CFLAGS)

clean:
	rm -f $(BUILDDIR)* $(OBJ_DIR)*.o
