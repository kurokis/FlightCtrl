#ifndef MOTORS_H_
#define MOTORS_H_


#include <inttypes.h>


enum BLCError
{
  BLC_ERROR_NONE = 0,
  BLC_ERROR_MISSING_MOTOR = 1<<0,
  BLC_ERROR_EXTRA_MOTOR = 1<<1,
  BLC_ERROR_INCONSISTENT_SETTINGS = 1<<2,
};


// =============================================================================
// Accessors

enum BLCError blc_error_bitfield(void);


// =============================================================================
// Public functions:

void DetectMotors(void);

// -----------------------------------------------------------------------------
void SetMotorSetpoint(uint8_t address, uint16_t setpoint);

// -----------------------------------------------------------------------------
void TxMotorSetpoints(void);


#endif  // MOTORS_H_
