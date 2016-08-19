/**\file */
#ifndef SLIC_DECLARATIONS_Poisson_H
#define SLIC_DECLARATIONS_Poisson_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define Poisson_PCIE_ALIGNMENT (16)
#define Poisson_L (32)
#define Poisson_M (32)
#define Poisson_N (32)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] ticks_PoissonKernel The number of ticks for which kernel "PoissonKernel" will run.
 * \param [in] inscalar_PoissonKernel_dh Input scalar parameter "PoissonKernel.dh".
 * \param [in] instream_poissonIn Stream "poissonIn".
 * \param [in] instream_size_poissonIn The size of the stream instream_poissonIn in bytes.
 * \param [out] outstream_poissonOut Stream "poissonOut".
 * \param [in] outstream_size_poissonOut The size of the stream outstream_poissonOut in bytes.
 * \param [in] inmem_PoissonKernel_twiddles Mapped ROM inmem_PoissonKernel_twiddles, should be of size (64 * sizeof(double)).
 */
void Poisson(
	uint64_t ticks_PoissonKernel,
	double inscalar_PoissonKernel_dh,
	const void *instream_poissonIn,
	size_t instream_size_poissonIn,
	void *outstream_poissonOut,
	size_t outstream_size_poissonOut,
	const double *inmem_PoissonKernel_twiddles);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] ticks_PoissonKernel The number of ticks for which kernel "PoissonKernel" will run.
 * \param [in] inscalar_PoissonKernel_dh Input scalar parameter "PoissonKernel.dh".
 * \param [in] instream_poissonIn Stream "poissonIn".
 * \param [in] instream_size_poissonIn The size of the stream instream_poissonIn in bytes.
 * \param [out] outstream_poissonOut Stream "poissonOut".
 * \param [in] outstream_size_poissonOut The size of the stream outstream_poissonOut in bytes.
 * \param [in] inmem_PoissonKernel_twiddles Mapped ROM inmem_PoissonKernel_twiddles, should be of size (64 * sizeof(double)).
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poisson_nonblock(
	uint64_t ticks_PoissonKernel,
	double inscalar_PoissonKernel_dh,
	const void *instream_poissonIn,
	size_t instream_size_poissonIn,
	void *outstream_poissonOut,
	size_t outstream_size_poissonOut,
	const double *inmem_PoissonKernel_twiddles);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t ticks_PoissonKernel; /**<  [in] The number of ticks for which kernel "PoissonKernel" will run. */
	double inscalar_PoissonKernel_dh; /**<  [in] Input scalar parameter "PoissonKernel.dh". */
	const void *instream_poissonIn; /**<  [in] Stream "poissonIn". */
	size_t instream_size_poissonIn; /**<  [in] The size of the stream instream_poissonIn in bytes. */
	void *outstream_poissonOut; /**<  [out] Stream "poissonOut". */
	size_t outstream_size_poissonOut; /**<  [in] The size of the stream outstream_poissonOut in bytes. */
	const double *inmem_PoissonKernel_twiddles; /**<  [in] Mapped ROM inmem_PoissonKernel_twiddles, should be of size (64 * sizeof(double)). */
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

