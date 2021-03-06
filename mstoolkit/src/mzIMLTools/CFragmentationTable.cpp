/*
Copyright 2017, Michael R. Hoopmann, Institute for Systems Biology
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "CFragmentationTable.h"

using namespace std;

CFragmentationTable::CFragmentationTable(){
  CMeasure m;
  measure = new vector<CMeasure>;
  measure->push_back(m);
}

CFragmentationTable::CFragmentationTable(const CFragmentationTable& c){
  measure = new vector<CMeasure>;
  for(size_t i=0;i<c.measure->size();i++) measure->push_back(c.measure->at(i));
}

CFragmentationTable::~CFragmentationTable(){
  delete measure;
}

CFragmentationTable& CFragmentationTable::operator=(const CFragmentationTable& c){
  if (this != &c){
    delete measure;
    measure = new vector<CMeasure>;
    for (size_t i = 0; i<c.measure->size(); i++) measure->push_back(c.measure->at(i));
  }
  return *this;
}

void CFragmentationTable::writeOut(FILE* f, int tabs){

}

