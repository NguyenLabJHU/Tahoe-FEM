/* $Id: PerTabEntryT.cpp,v 1.1 2002-03-13 02:25:55 jzimmer Exp $ */

#include "PerTabEntryT.h"

#include "dArrayT.h"
#include "StringT.h"

PerTabEntryT::PerTabEntryT() {

        atomicmass = 0;
}

PerTabEntryT::~PerTabEntryT() {
}

void PerTabEntryT::SetName(const char * s) {
	name = s;
}
void PerTabEntryT::SetSymbol(const char * s) {
	symbol = s;
}	

void PerTabEntryT::SetLattParams(double a) {
        latt_params[0] = latt_params[1] = latt_params[2] = a;
}

void PerTabEntryT::SetLattParams(double a, double c) {
        latt_params[0] = latt_params[1] = a;
	latt_params[2] = c;
}

void PerTabEntryT::SetLattParams(double a, double b, double c) {
	latt_params[0] = a;
	latt_params[1] = b;
	latt_params[2] = c;
}
 
void PerTabEntryT::SetMass(double m) {
	atomicmass = m;
}

void PerTabEntryT::SetLattType(const char * s) {
	latticetype = s;
} 

const StringT& GetName() {
	return name;
}

const StringT& GetSymbol() {
	return symbol;
}

const dArrayT& GetLattParams() {
	return latt_params;
}

const double& GetMass() {
	return atomicmass;
}

const StringT& GetLattType() {
	return latticetype;
}

