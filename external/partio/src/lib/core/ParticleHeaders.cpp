/*
PARTIO SOFTWARE
Copyright 2010 Disney Enterprises, Inc. All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.

* The names "Disney", "Walt Disney Pictures", "Walt Disney Animation
Studios" or the names of its contributors may NOT be used to
endorse or promote products derived from this software without
specific prior written permission from Walt Disney Pictures.

Disclaimer: THIS SOFTWARE IS PROVIDED BY WALT DISNEY PICTURES AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, NONINFRINGEMENT AND TITLE ARE DISCLAIMED.
IN NO EVENT SHALL WALT DISNEY PICTURES, THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND BASED ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
*/
#include "ParticleHeaders.h"
#include <map>
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace Partio;

ParticleHeaders::
ParticleHeaders()
    :particleCount(0)
{
}

ParticleHeaders::
~ParticleHeaders()
{
}

void ParticleHeaders::
release()
{
    delete this;
}

int ParticleHeaders::
nuParticles() const
{
    return particleCount;
}

int ParticleHeaders::
numAttributes() const
{
    return static_cast<int>(attributes.size());
}

int ParticleHeaders::
numFixedAttributes() const
{
    return static_cast<int>(fixedAttributes.size());
}

bool ParticleHeaders::
attributeInfo(const int attributeIndex,ParticleAttribute& attribute) const
{
    if(attributeIndex<0 || attributeIndex>=(int)attributes.size()) return false;
    attribute=attributes[attributeIndex];
    return true;
}

bool ParticleHeaders::
fixedAttributeInfo(const int attributeIndex,FixedAttribute& attribute) const
{
    if(attributeIndex<0 || attributeIndex>=(int)fixedAttributes.size()) return false;
    attribute=fixedAttributes[attributeIndex];
    return true;
}

bool ParticleHeaders::
attributeInfo(const char* attributeName,ParticleAttribute& attribute) const
{
    std::map<std::string,int>::const_iterator it=nameToAttribute.find(attributeName);
    if(it!=nameToAttribute.end()){
        attribute=attributes[it->second];
        return true;
    }
    return false;
}

bool ParticleHeaders::
fixedAttributeInfo(const char* attributeName,FixedAttribute& attribute) const
{
    std::map<std::string,int>::const_iterator it=nameToFixedAttribute.find(attributeName);
    if(it!=nameToFixedAttribute.end()){
        attribute=fixedAttributes[it->second];
        return true;
    }
    return false;
}

void ParticleHeaders::
sort()
{
    assert(false);
}


int ParticleHeaders::
registerIndexedStr(const ParticleAttribute&, const char*)
{
    assert(false);
    return -1;
}

int ParticleHeaders::
registerFixedIndexedStr(const FixedAttribute&, const char* str)
{
    assert(false);
    return -1;
}

int ParticleHeaders::
lookupIndexedStr(const ParticleAttribute&, const char*) const
{
    assert(false);
    return -1;
}

int ParticleHeaders::
lookupFixedIndexedStr(const FixedAttribute&, const char* str) const
{
    assert(false);
    return -1;
}

const std::vector<std::string>& ParticleHeaders::
indexedStrs(const ParticleAttribute&) const
{
    static std::vector<std::string> dummy;
    assert(false);
    return dummy;
}

const std::vector<std::string>& ParticleHeaders::
fixedIndexedStrs(const FixedAttribute& attr) const
{
    static std::vector<std::string> dummy;
    assert(false);
    return dummy;
}

void ParticleHeaders::
findPoints(const float[3],const float[3],std::vector<ParticleIndex>&) const
{
    assert(false);
}

float ParticleHeaders::
findNPoints(const float[3],const int,const float,std::vector<ParticleIndex>&,std::vector<float>&) const
{
    assert(false);
    return 0;
}

int ParticleHeaders::
findNPoints(const float[3],int,const float, ParticleIndex *,
    float *, float *) const
{
    assert(false);
    return 0;
}

ParticleAttribute ParticleHeaders::
addAttribute(const char* attribute,ParticleAttributeType type,const int count)
{
    // TODO: check if attribute already exists and if so what data type
    ParticleAttribute attr;
    attr.name=attribute;
    attr.type=type;
    attr.attributeIndex=static_cast<int>(attributes.size()); //  all arrays separate so we don't use this here!
    attr.count=count;
    attributes.push_back(attr);
    nameToAttribute[attribute]=static_cast<int>(attributes.size()-1);
    return attr;
}

FixedAttribute ParticleHeaders::
addFixedAttribute(const char* attribute,ParticleAttributeType type,const int count)
{
    // TODO: check if attribute already exists and if so what data type
    FixedAttribute attr;
    attr.name=attribute;
    attr.type=type;
    attr.attributeIndex=fixedAttributes.size(); //  all arrays separate so we don't use this here!
    attr.count=count;
    fixedAttributes.push_back(attr);
    nameToFixedAttribute[attribute]=fixedAttributes.size()-1;
    return attr;
}

ParticleIndex ParticleHeaders::
addParticle()
{
    ParticleIndex index=particleCount;
    particleCount++;
    return index;
}

ParticlesDataMutable::iterator ParticleHeaders::
addParticles(const int countToAdd)
{
    particleCount+=countToAdd;
    return iterator();
}

void* ParticleHeaders::
dataInternal(const ParticleAttribute&,const ParticleIndex) const
{
    assert(false);
    return 0;
}

void* ParticleHeaders::
fixedDataInternal(const FixedAttribute& attribute) const
{
    assert(false);
    return 0;
}

void ParticleHeaders::
dataInternalMultiple(const ParticleAttribute&,const int,
    const ParticleIndex*,const bool,char*) const
{
    assert(false);
}

void ParticleHeaders::
dataAsFloat(const ParticleAttribute&,const int,
    const ParticleIndex*,const bool,float*) const
{
    assert(false);
}


void ParticleHeaders::
setIndexedStr(const ParticleAttribute& attribute,int indexedStringToken,const char* str){
    assert(false);
}

void ParticleHeaders::
setFixedIndexedStr(const FixedAttribute& attribute,int indexedStringToken,const char* str){
    assert(false);
}
