#include<stdlib.h>
#include<string.h>
#include<stdio.h>

/**
*  inserts a character string after the filename, before the extension
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void name_ext(char *filename, char *extension)
{
	char ext[100];

	if (strstr(filename, ".su") != NULL) {
		sprintf(ext,"%s.su", extension);
		strcpy(strstr(filename, ".su"), ext);
	}
	else if (strstr(filename, ".segy") != NULL) {
		sprintf(ext,"%s.segy", extension);
		strcpy(strstr(filename, ".segy"), ext);
	}
	else if (strstr(filename, ".mat") != NULL) {
		sprintf(ext,"%s.mat", extension);
		strcpy(strstr(filename, ".mat"), ext);
	}
	else if (strstr(filename, ".hdf") != NULL) {
		sprintf(ext,"%s.hdf", extension);
		strcpy(strstr(filename, ".hdf"), ext);
	}
	else if (strrchr(filename, '.') != NULL) {
		sprintf(ext,"%s.su", extension);
		strcpy(strrchr(filename, '.'), ext);
	}
	else {
		sprintf(ext,"%s.su", extension);
		strcat(filename, ext);
	}

	return;
}
