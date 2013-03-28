/*
 * ReactionEquation.cpp
 *
 *  Created on: 9 Nov 2012
 *      Author: robinsonm
 */

#include "ReactionEquation.h"

namespace Tyche {
ReactionComponent operator*(const int mult, Species& s) {
	return ReactionComponent(mult,s,0);
}

ReactionSide operator+(const ReactionComponent& arg1, const ReactionComponent& arg2) {
	return ReactionSide(arg1,arg2);
}

ReactionSide operator+(Species& arg1, const ReactionComponent& arg2) {
	return ReactionSide(ReactionComponent(1,arg1,0),arg2);
}

ReactionSide operator+(const ReactionComponent& arg1, Species& arg2) {
	return ReactionSide(arg1,ReactionComponent(1,arg2,0));
}

ReactionSide operator+(Species& arg1, Species& arg2) {
	return ReactionSide(ReactionComponent(1,arg1,0),ReactionComponent(1,arg2,0));
}

ReactionSide& operator+(ReactionSide& side, const ReactionComponent& comp) {
	side.push_back(comp);
	return side;
}

ReactionSide& operator+(ReactionSide& side, Species& s) {
	side.push_back(ReactionComponent(1,s,0));
	return side;
}

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

ReactionEquation operator>>(const ReactionSide& lhs, Species& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs,0));
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
ReactionEquation operator>>(const ReactionComponent& lhs, Species& rhs) {
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs,0));
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(const ReactionComponent& lhs, const int rhs) {
	ASSERT(rhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide(lhs);
	ReactionSide* newrhs = new ReactionSide();
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(Species& lhs, const ReactionSide& rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs,0));
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}

ReactionEquation operator>>(Species& lhs, const ReactionComponent& rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs,0));
	ReactionSide* newrhs = new ReactionSide(rhs);
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(Species& lhs, Species& rhs) {
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs,0));
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs,0));
	return ReactionEquation(*newlhs,*newrhs);
}
ReactionEquation operator>>(Species& lhs, const int rhs) {
	ASSERT(rhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide(ReactionComponent(1,lhs,0));
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
ReactionEquation operator>>(const int lhs, Species& rhs) {
	ASSERT(lhs==0,"null species is always 0");
	ReactionSide* newlhs = new ReactionSide();
	ReactionSide* newrhs = new ReactionSide(ReactionComponent(1,rhs,0));
	return ReactionEquation(*newlhs,*newrhs);
}


ReactionEquation unimolecular_reaction(const int i1, Species& s1, const int i2, Species& s2) {
	ReactionSide* lhs = new ReactionSide();
	ReactionSide* rhs = new ReactionSide();
	if (s1.id != null_species.id) {
		lhs->push_back(ReactionComponent(i1,s1,0.0));
	}
	if (s2.id != null_species.id) {
		rhs->push_back(ReactionComponent(i2,s2,0.0));
	}
	return ReactionEquation(*lhs,*rhs);
}

ReactionEquation bimolecular_reaction(const int i1, Species& s1, const int i2, Species& s2, const int i3, Species& s3) {
	ReactionSide* lhs = new ReactionSide();
	ReactionSide* rhs = new ReactionSide();
	if (s1.id == s2.id) {
		lhs->push_back(ReactionComponent(i1+i2,s1,0.0));
	} else {
		lhs->push_back(ReactionComponent(i1,s1,0.0));
		lhs->push_back(ReactionComponent(i2,s2,0.0));
	}
	if (s3.id != null_species.id) {
		rhs->push_back(ReactionComponent(i3,s3,0.0));
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
