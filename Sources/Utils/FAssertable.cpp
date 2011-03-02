#include "FAssertable.hpp"
#include "FAbstractApplication.hpp"

/**
* Current application
*/
FAbstractApplication* FAssertable::CurrentApp(0);

// Simply quit the current app
void FAssertable::exitApplication(const int inExitCode) const{
	if(CurrentApp) CurrentApp->abort(inExitCode);
}




