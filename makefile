CC        =g++
CFLAGS    =-c -std=gnu++11 -Wall -O3 -I/usr/local/Cellar/gsl/2.3/include/gsl
LDFLAGS   =-lstdc++ -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas -lm
INCLUDE   =-I./include/  -I/usr/local/Cellar/gsl/1.16/include/gsl
OBJDIR    =obj/

# Backslash for linebreak. Wohoo!! Yay! (and comment is # obv. Beware, it affects the *whole* line.)
OBJLIST   = Point.o Farm.o Shipment_kernel.o \
		 	Region.o County.o State.o Grid_cell.o Grid_checker.o \
			Grid_manager.o Status_manager.o	shared_functions.o \
			File_manager.o main.o Control_manager.o Shipment_manager.o \
			Local_spread.o Control_resource.o USAMM_parameters.o Population_manager.o

OBJECTS   = $(addprefix $(OBJDIR), $(OBJLIST) )

all:USDOSv2_new

USDOSv2_new: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OBJECTS): ./$(OBJDIR)%.o: src/%.cpp
	$(CC) $(CFLAGS) $? -o $@ $(INCLUDE)

clean:
	rm obj/*.o
