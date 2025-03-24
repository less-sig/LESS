#include <stdint.h>

#include "stm32f4xx_hal.h"

#if defined(STM32F469xx)
#define LED_PIN GPIO_PIN_4
#elif defined(STM32F412Zx)
#define LED_PIN GPIO_PIN_3
#else
#define LED_PIN GPIO_PIN_12
#endif

#if defined(STM32F412Zx)
#define LED_GPIO_PORT GPIOE
#define LED_GPIO_CLK_ENABLE() __HAL_RCC_GPIOE_CLK_ENABLE()
#else
#define LED_GPIO_PORT GPIOD
#define LED_GPIO_CLK_ENABLE() __HAL_RCC_GPIOD_CLK_ENABLE()
#endif

void utils_init(void);
void LED_toggle(void);
void stack_usage_scan_init(uint32_t max_stack_size);
uint32_t stack_usage_scan(void);

#define bm_decls uint32_t bm_s, bm_e
#define bm_start() bm_s = DWT->CYCCNT
#define bm_end() bm_e = DWT->CYCCNT
#define bm_result() (bm_e - bm_s)

static inline
uint32_t get_cycles(void) {
	return DWT->CYCCNT;
}

extern volatile uint32_t stack_usage_sample_current, stack_usage_sample_max;

#define stack_usage_sample_get_last() stack_usage_sample_current
#define stack_usage_sample_get_max() stack_usage_sample_max

#define FIX_RAM_FUNC(func) typeof(func) *func##_fixed = (typeof(func) *)((uint32_t)func & 0xDFFFFFFF)
