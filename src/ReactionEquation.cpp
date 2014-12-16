/*
 * ReactionEquation.cpp
 *
 *  Created on: 9 Nov 2012
 *      Author: robinsonm
 */

#include "ReactionEquation.h"

namespace Tyche {
ReactionComponent operator*(const int mult, SpeciesType s) {
	return ReactionComponent(mult,s.s,s.type);
}

ReactionSide operator+(const ReactionComponent& arg1, const ReactionComponent& arg2) {
	return ReactionSide(arg1,arg2);
}

ReactionSide operator+(SpeciesType arg1, const ReactionComponent& arg2) {
	return ReactionSide(ReactionComponent(1,arg1.s,arg1.type),arg2);
}

ReactionSide operator+(const ReactionComponent& arg1, SpeciesType arg2) {
	return ReactionSide(arg1,ReactionComponent(1,arg2.s,arg2.type));
}

ReactionSide operator+(SpeciesType arg1, SpeciesType arg2) {
	return ReactionSide(ReactionComponent(1,arg1.s,arg1.type),ReactionComponent(1,arg2.s,arg2.type));
}

//ReactionSide& operator+(ReactionSide& side, const ReactionComponent& comp) {
//	side.push_back(comp);
//	return side;
//}

ReactionSide operator+(const ReactionSide side, const ReactionComponent& comp) {
	ReactionSide result(side);
	result.push_back(comp);
	return result;
}

ReactionSide operator+(const ReactionSide side, SpeciesType s) {
	ReactionSide result(side);
	result.push_back(ReactionComponent(1,s.s,s.type));
	return result;
}

//ReactionSide& operator+(ReactionSide& side, SpeciesType s) {
//	side.push_back(ReactionComponent(1,s,0));
//	return side;
//}

ReactionEquation operator>>(const ReactionSide& lhs, const ReactionSide& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const ReactionSide& lhs, const ReactionComponent& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const ReactionSide& lhs, SpeciesType rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs.s,rhs.type));
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const ReactionSide& lhs, const int rhs) {
	ASSERT(rhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide();
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const ReactionComponent& lhs, const ReactionSide& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const ReactionComponent& lhs, const ReactionComponent& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(const ReactionComponent& lhs, SpeciesType rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs.s,rhs.type));
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(const ReactionComponent& lhs, const int rhs) {
	ASSERT(rhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide();
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(SpeciesType lhs, const ReactionSide& rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs.s,lhs.type));
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(SpeciesType lhs, const ReactionComponent& rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs.s,lhs.type));
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(SpeciesType lhs, SpeciesType rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs.s,lhs.type));
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs.s,rhs.type));
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(SpeciesType lhs, const int rhs) {
	ASSERT(rhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs.s,lhs.type));
	ReactionSide* newrhs = new ReactionSide();
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const int lhs, const ReactionSide& rhs) {
	ASSERT(lhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide();
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(const int lhs, const ReactionComponent& rhs) {
	ASSERT(lhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide();
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(const int lhs, SpeciesType rhs) {
	ASSERT(lhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide();
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs.s,rhs.type));
	return ReactionEquation(*newlhs,*newrhs);
}


ReactionEquation unimolecular_reaction(const int i1, SpeciesType s1, const int i2, SpeciesType s2) {
	ReactionSide* lhs = new ReactionSide();
	ReactionSide* rhs = new ReactionSide();
	if (s1.s->id != null_species.id) {
		lhs->push_back(ReactionComponent(i1,s1.s,s1.type));
	}
	if (s2.s->id != null_species.id) {
		rhs->push_back(ReactionComponent(i2,s2.s,s2.type));
	}
	return ReactionEquation(*lhs,*rhs);
}

ReactionEquation bimolecular_reaction(const int i1, SpeciesType s1, const int i2, SpeciesType s2, const int i3, SpeciesType s3) {
	ReactionSide* lhs = new ReactionSide();
	ReactionSide* rhs = new ReactionSide();
	if (s1.s->id == s2.s->id) {
		lhs->push_back(ReactionComponent(i1+i2,s1.s,s1.type));
	} else {
		lhs->push_back(ReactionComponent(i1,s1.s,s1.type));
		lhs->push_back(ReactionComponent(i2,s2.s,s2.type));
	}
	if (s3.s->id != null_species.id) {
		rhs->push_back(ReactionComponent(i3,s3.s,s3.type));
	}
	return ReactionEquation(*lhs,*rhs);
}

std::ostream& operator<< (std::ostream& out, const ReactionSide &side) {
	const int n = side.size();
	for (int i = 0; i < n; ++i) {
		out << side[i].multiplier << "("<<side[i].species->id<<")";
		if (i != n-1) {
			out << " + ";
		}
	}
	return out;
}
std::ostream& operator<< (std::ostream& out, const ReactionEquation &eq) {
	return out << eq.lhs << " >> " << eq.rhs;
}
}
