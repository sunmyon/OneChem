#ifndef STANDALONE_HOST_H
#define STANDALONE_HOST_H

/* Process termination is owned by the executable, never by portable core. */
void endrun(void);

/* Compatibility values required only by the pre-portable chemistry modules. */
#if !defined(CHEMISTRY_PORTABLE) && !defined(CHEMISTRY_OLD)
extern const int ThisTask;
extern const int NTask;
#endif

#endif
