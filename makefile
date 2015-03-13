TARGET := $(notdir $(shell pwd))

MCU    := atmega1284p
F_CPU  := 20000000

CCFLAGS    = -std=gnu99 -Wstrict-prototypes
LSTFLAGS   = -Wa,-adhlns=$(LST)
# Temporarily removed -Werror until gcc 4.8.3 comes to Ubuntu
LDFLAGS    = -flto -Ofast -fwhole-program -pedantic -Wall -Wextra -Wundef \
             -fshort-enums -ffreestanding -ffunction-sections -fdata-sections \
             -Wl,--relax,--gc-sections
ALLFLAGS   = -mmcu=$(MCU) -DF_CPU="$(F_CPU)UL"
DUDEFLAGS  = -c avrisp2 -p $(MCU)

CC     := avr-gcc
CP     := avr-objcopy
DUMP   := avr-objdump
DUDE   := avrdude

# If the environment variable DEV_BUILD_PATH is set, then the build files will
# be placed there in a named sub-folder, otherwise a build directory will be
# created in the current directory
ifneq ($(DEV_BUILD_PATH),)
  BUILD_PATH := $(DEV_BUILD_PATH)/build/$(TARGET)
else
  BUILD_PATH := build
endif

SOURCES  := $(wildcard *.c)
SOURCES  += $(wildcard *.S)
ASSEMBLY := $(addsuffix .lst, $(addprefix $(BUILD_PATH)/, $(SOURCES)))

ELF    := $(BUILD_PATH)/$(TARGET).elf
HEX    := $(BUILD_PATH)/$(TARGET).hex
LST    := $(BUILD_PATH)/$(TARGET).lst

# Rules to make the assembly listings
$(BUILD_PATH)/%.c.lst: %.c
	$(CC) -c $(LDFLAGS) $(CCFLAGS) $(ALLFLAGS) -Wa,-adhlns=$@ -o /dev/null $<

$(BUILD_PATH)/%.S.lst: %.S
	$(CC) -c $(LDFLAGS) $(CCFLAGS) $(ALLFLAGS) -Wa,-adhlns=$@ -o /dev/null $<

# Declare targets that are not files
.PHONY: program clean

all: $(HEX) $(LST)

$(HEX): $(ELF)
	$(CP) -O ihex -R .eeprom $< $@

# Target to make assembly listing of link-time full-program optimized output.
$(LST): $(ELF)
	$(DUMP) -d $(ELF) > $(LST)

# Target to build the .elf file
# NOTE: -lm includes the math library (libm.a)
$(ELF): $(SOURCES) $(BUILD_PATH)
	$(CC) $(LDFLAGS) $(CCFLAGS) $(ALLFLAGS) $(LSTFLAGS) -o $@ $(SOURCES) -lm

# Target to program the microprocessor flash only
program: $(HEX)
	$(DUDE) $(DUDEFLAGS) -U flash:w:$(HEX):i

# Target to make assembly listings.
# WARNING!!!: Because this makefile employs link-time optimization, the final
# program may be different from what appears in these files!!!
# Listings are first cleared here since there are no dependency checks.
assembly: clean_assembly $(BUILD_PATH) $(ASSEMBLY)

# Target to clean up the directory (leaving only source)
clean: $(BUILD_PATH)
	rm -f $(HEX) $(ELF) $(LST) $(ASSEMBLY)
	rmdir $(BUILD_PATH)

clean_assembly:
	rm -f $(ASSEMBLY)

$(BUILD_PATH):
	mkdir -p $(BUILD_PATH)