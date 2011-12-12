

benchmark-arm.elf: benchmark.c interface-arm.c ../fix16_base.c
	arm-none-eabi-gcc -mcpu=cortex-m3 -mthumb -T generic-m-hosted.ld \
		-Wall -O2 \
		-o $@ -I .. $^

run-benchmark-arm: benchmark-arm.elf
	qemu-system-arm -cpu cortex-m3 -icount 0 -device armv7m_nvic \
		-nographic -monitor null -serial null \
		-semihosting -kernel $<

benchmark-avr.elf: benchmark.c interface-avr.c ../fix16_base.c
	avr-gcc -Wall -mmcu=atmega128 \
		-Wall -O2 -DFIXMATH_NO_64BIT \
		-o $@ -I .. $^

run-benchmark-avr: benchmark-avr.elf
	# Note: this needs simulavrxx 1.0rc0 or newer
	simulavr -d atmega128 -f $< -W 0x20,- -T exit
