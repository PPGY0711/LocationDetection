#include "errormsg.h"

void errorReport(char* errormsg, ...)
{
	va_list args;
	fprintf(stderr, "ERROR: ");
	va_start(args, errormsg);
	vfprintf(stderr, errormsg, args);
	va_end(args);
	fprintf(stderr, "\n");
	return;
}
