/* $Id: PerTabEntryT.cpp,v 1.2 2002-06-03 22:20:43 jzimmer Exp $ */

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

const StringT& PerTabEntryT::GetName() {
	return name;
}

const StringT& PerTabEntryT::GetSymbol() {
	return symbol;
}

const dArrayT& PerTabEntryT::GetLattParams() {
	return latt_params[1];
}

const double& PerTabEntryT::GetMass() {
	return atomicmass;
}

const StringT& PerTabEntryT::GetLattType() {
	return latticetype;
}

