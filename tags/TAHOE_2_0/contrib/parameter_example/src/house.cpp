/* $Id: house.cpp,v 1.4 2003-11-21 22:48:26 paklein Exp $ */
#include "house.h"
#include "lawn.h"
#include "AutoArrayT.h"
#include "StringT.h"

/* rooms */
#include "closet.h"
#include "bedroom.h"

/* basements */
#include "crawl_space.h"
#include "storm_shelter.h"

house::house(void):
	ParameterInterfaceT("house"),
	style_(undefined),
	zipcode_(-1),
	basement_(NULL),
	lawn_(NULL)
{

}

house::~house(void) {
	for (int i = 0; i < rooms_.Length(); i++)
		delete rooms_[i];
	delete lawn_;
	delete basement_;
}

void house::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	ParameterT the_style(ParameterT::Enumeration, "style");
	the_style.AddEnumeration("colonial", colonial);
	the_style.AddEnumeration("ranch", ranch);
	the_style.AddEnumeration("split_level", split_level);
	list.AddParameter(the_style);

	list.AddParameter(zipcode_, "zipcode");	
}

void house::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "house::TakeParameterList";

	/* my parameters */
	list.GetParameter("style", enum2int<style>(style_));
	zipcode_ = list.GetParameter("zipcode");

	/* the roof */
	const ParameterListT* roof_params = list.List("roof");
	if (!roof_params) ExceptionT::GeneralFail(caller, "cannot find \"roof\"");
	roof_.TakeParameterList(*roof_params);

	/* the driveway */
	const ParameterListT* driveway_params = list.List("driveway");
	if (!driveway_params) ExceptionT::GeneralFail(caller, "cannot find \"driveway\"");
	driveway_.TakeParameterList(*driveway_params);
	
	/* 1st garage */
	const ParameterListT* garage_params_1 = list.List("garage", 1);
	if (!garage_params_1) ExceptionT::GeneralFail(caller, "cannot find 1st \"garage\"");
	garage1_.TakeParameterList(*garage_params_1);

	/* 2nd garage */
	const ParameterListT* garage_params_2 = list.List("garage", 2);
	if (!garage_params_2) ExceptionT::GeneralFail(caller, "cannot find 2nd \"garage\"");
	garage2_.TakeParameterList(*garage_params_2);

	/* (optional) lawn */
	const ParameterListT* lawn_params = list.List("lawn");
	if (lawn_params) {
		lawn_ = new lawn;
		lawn_->TakeParameterList(*lawn_params);
	}

	/* (optional) basement */
	const ParameterListT* basement_params = list.List("crawl_space");
	if (!basement_params) basement_params = list.List("storm_shelter");
	if (basement_params) {
		basement_ = New_basement(basement_params->Name(), true);
		basement_->TakeParameterList(*basement_params);
	}
	
	/* rooms */
	int num_closets = list.NumLists("closet");
	int num_bedrooms = list.NumLists("bedroom");
	int num_rooms = num_closets + num_bedrooms;
	rooms_.Dimension(num_rooms);
	num_rooms = 0;
	for (int i = 1; i <= num_closets; i++) {
		const ParameterListT* closet_params = list.List("closet", i);
		if (!closet_params) ExceptionT::GeneralFail(caller, "cannot find \"closet\" %d", i);
		room* new_room = New_room(closet_params->Name(), true);
		new_room->TakeParameterList(*closet_params);
		rooms_[num_rooms++] = new_room;
	}
	for (int i = 1; i <= num_bedrooms; i++) {
		const ParameterListT* br_params = list.List("bedroom", i);
		if (!br_params) ExceptionT::GeneralFail(caller, "cannot find \"bedroom\" %d", i);
		room* new_room = New_room(br_params->Name(), true);
		new_room->TakeParameterList(*br_params);
		rooms_[num_rooms++] = new_room;
	}
}

void house::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParameterInterfaceT::DefineSubs(sub_list);

	/* the roof */
	sub_list.AddSub("roof");

	/* the driveway */
	sub_list.AddSub("driveway");

	/* 2 garages */
	SubListDescriptionT garage("garage");
	sub_list.AddSub(garage);
	sub_list.AddSub(garage);

	/* the lawn */
	sub_list.AddSub("lawn", ParameterListT::ZeroOrOnce);

	/* the basement */
	sub_list.AddSub("basement", ParameterListT::ZeroOrOnce, true);

	/* the rooms */
	sub_list.AddSub("rooms", ParameterListT::OnePlus, true);
}

void house::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
	SubListT& sub_sub_list) const
{
	if (sub == "rooms")
	{
		order = ParameterListT::Choice;
		
		/* a closet */
		sub_sub_list.AddSub("closet");

		/* a bedroom */
		sub_sub_list.AddSub("bedroom");
	}
	else if (sub == "basement")
	{
		order = ParameterListT::Choice;
		
		/* a crawl space */
		sub_sub_list.AddSub("crawl_space");

		/* a bedroom */
		sub_sub_list.AddSub("storm_shelter");
	}
	else /* inherited */
		ParameterInterfaceT::DefineInlineSub(sub, order, sub_sub_list);
}

ParameterInterfaceT* house::NewSub(const StringT& list_name) const
{
	/* try rooms */
	room* new_room = New_room(list_name, false);
	if (new_room) return new_room;
	
	/* try basements */
	basement* new_basement = New_basement(list_name, false);
	if (new_basement) return new_basement;

	/* try others */
	if (list_name == "roof") {
		return new roof;
	}
	else if (list_name == "driveway") {
		return new driveway;
	}
	else if (list_name == "lawn") {
		return new lawn;
	}
	else if (list_name == "garage") {
		return new garage;
	}
	else /* inherited */
		return ParameterInterfaceT::NewSub(list_name);
}

/**********************************************************************
 * Private
 **********************************************************************/

room* house::New_room(const StringT& room_name, bool throw_on_fail) const
{
	if (room_name == "closet")
		return new closet;
	else if (room_name == "bedroom")
		return new bedroom;
	else if (throw_on_fail)
		ExceptionT::GeneralFail("house::New_room", "unrecognized room type \"%s\"", room_name.Pointer());
	return NULL;
}

basement* house::New_basement(const StringT& basement_name, bool throw_on_fail) const
{
	if (basement_name == "crawl_space")
		return new crawl_space;
	else if (basement_name == "storm_shelter")
		return new storm_shelter;
	else if (throw_on_fail)
		ExceptionT::GeneralFail("house::New_basement", "unrecognized basement type \"%s\"", basement_name.Pointer());
	
	return NULL;
}
