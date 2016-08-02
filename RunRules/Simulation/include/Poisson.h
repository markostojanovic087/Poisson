/**\file */
#ifndef SLIC_DECLARATIONS_Poisson_H
#define SLIC_DECLARATIONS_Poisson_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define Poisson_PCIE_ALIGNMENT (16)
#define Poisson_N (64)


/*----------------------------------------------------------------------------*/
/*-------------------------- Interface MyInterface ---------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'MyInterface'.
 * 
 * \param [in] param_N Interface Parameter "N".
 * \param [in] instream_fft_in The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes.
 * \param [out] outstream_fft_out The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes.
 */
void Poisson_MyInterface(
	uint64_t param_N,
	const float *instream_fft_in,
	float *outstream_fft_out);

/**
 * \brief Basic static non-blocking function for the interface 'MyInterface'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_N Interface Parameter "N".
 * \param [in] instream_fft_in The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes.
 * \param [out] outstream_fft_out The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poisson_MyInterface_nonblock(
	uint64_t param_N,
	const float *instream_fft_in,
	float *outstream_fft_out);

/**
 * \brief Advanced static interface, structure for the engine interface 'MyInterface'
 * 
 */
typedef struct { 
	uint64_t param_N; /**<  [in] Interface Parameter "N". */
	const float *instream_fft_in; /**<  [in] The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes. */
	float *outstream_fft_out; /**<  [out] The stream should be of size ((((param_N * 2) * param_N) * param_N) * 4) bytes. */
} Poisson_MyInterface_actions_t;

/**
 * \brief Advanced static function for the interface 'MyInterface'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Poisson_MyInterface_run(
	max_engine_t *engine,
	Poisson_MyInterface_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'MyInterface'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_MyInterface_run_nonblock(
	max_engine_t *engine,
	Poisson_MyInterface_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'MyInterface'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Poisson_MyInterface_run_group(max_group_t *group, Poisson_MyInterface_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'MyInterface'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_MyInterface_run_group_nonblock(max_group_t *group, Poisson_MyInterface_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'MyInterface'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Poisson_MyInterface_run_array(max_engarray_t *engarray, Poisson_MyInterface_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'MyInterface'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_MyInterface_run_array_nonblock(max_engarray_t *engarray, Poisson_MyInterface_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Poisson_MyInterface_convert(max_file_t *maxfile, Poisson_MyInterface_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_PoissonKernel The number of ticks for which kernel "PoissonKernel" will run.
 * \param [in] instream_fft_in Stream "fft_in".
 * \param [in] instream_size_fft_in The size of the stream instream_fft_in in bytes.
 * \param [out] outstream_fft_out Stream "fft_out".
 * \param [in] outstream_size_fft_out The size of the stream outstream_fft_out in bytes.
 */
void Poisson(
	uint64_t ticks_PoissonKernel,
	const void *instream_fft_in,
	size_t instream_size_fft_in,
	void *outstream_fft_out,
	size_t outstream_size_fft_out);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] ticks_PoissonKernel The number of ticks for which kernel "PoissonKernel" will run.
 * \param [in] instream_fft_in Stream "fft_in".
 * \param [in] instream_size_fft_in The size of the stream instream_fft_in in bytes.
 * \param [out] outstream_fft_out Stream "fft_out".
 * \param [in] outstream_size_fft_out The size of the stream outstream_fft_out in bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poisson_nonblock(
	uint64_t ticks_PoissonKernel,
	const void *instream_fft_in,
	size_t instream_size_fft_in,
	void *outstream_fft_out,
	size_t outstream_size_fft_out);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_PoissonKernel; /**<  [in] The number of ticks for which kernel "PoissonKernel" will run. */
	const void *instream_fft_in; /**<  [in] Stream "fft_in". */
	size_t instream_size_fft_in; /**<  [in] The size of the stream instream_fft_in in bytes. */
	void *outstream_fft_out; /**<  [out] Stream "fft_out". */
	size_t outstream_size_fft_out; /**<  [in] The size of the stream outstream_fft_out in bytes. */
} Poisson_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Poisson_run(
	max_engine_t *engine,
	Poisson_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_run_nonblock(
	max_engine_t *engine,
	Poisson_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Poisson_run_group(max_group_t *group, Poisson_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_run_group_nonblock(max_group_t *group, Poisson_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Poisson_run_array(max_engarray_t *engarray, Poisson_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poisson_run_array_nonblock(max_engarray_t *engarray, Poisson_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Poisson_convert(max_file_t *maxfile, Poisson_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Poisson_init(void);

/* Error handling functions */
int Poisson_has_errors(void);
const char* Poisson_get_errors(void);
void Poisson_clear_errors(void);
/* Free statically allocated maxfile data */
void Poisson_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int Poisson_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int Poisson_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Poisson_H */

