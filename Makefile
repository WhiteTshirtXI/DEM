PACKAGE = PD

MEMBERS = main initialize calculate_contact_force contact_detect erase_contact_data get_cell_neighbors move_particle update_velocity sort_particles wall_position wall_velocity handle_walls get_material_properties

ifeq ($(DEBUG),yes)
CFLAGS = -g -pg 
LDFLAGS = -g -pg -lm
else
CFLAGS  = -D_GNU_SOURCE -O3 -ffast-math -fno-trapping-math -funroll-loops -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -MD
LDFLAGS = -Wl,-stack_size,0xf000000 -lm
endif

.SUFFIXES: .d

OBJS = $(patsubst %,%.o,$(MEMBERS))
DEPS = $(patsubst %,%.d,$(MEMBERS))
SOURCES = $(patsubst %,%.c,$(MEMBERS))

all: $(PACKAGE)

$(PACKAGE): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(PACKAGE) $(OBJS) $(DEPS)

rclean:
	rm -f $(PACKAGE) $(OBJS) $(DEPS) *.bck *.bak *.BAK *~

rrclean:
	rm -f $(PACKAGE) $(OBJS) $(DEPS) jam_times *.bck *.bak *.BAK *~ file* fluid*

-include $(DEPS)
