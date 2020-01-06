/*
Name: Utkarsh Srivastava
CWID:10832610
*/
//-------------------------------------Header Files-------------------------------------------

#include <stdio.h>
#include <ctype.h>
#define string_length strlen
#include <vector>
#include <string.h>
#include <algorithm>
#include "check_true_false.h"
#include <iostream>

using namespace std;

logical_expression::logical_expression()
{
  symbol[0] = 0;
  connective[0] = 0;
}


logical_expression::~logical_expression()
{
  unsigned long counter;
  for (counter = 0; counter < subexpressions.size(); counter++)
  {
    delete subexpressions[counter];
  }
}


void print_expression(logical_expression * expression, char * separator)
{
  if (expression == 0)
  {
    printf("\nINVALID\n");
  }
  else if (strcmp(expression->symbol, "") != 0)
  {
    printf("%s", expression->symbol);
  }
  else
  {
    printf("(%s",  expression->connective);
    unsigned long counter;
    for (counter = 0; counter < expression->subexpressions.size(); counter++)
    {
      printf(" ");
      print_expression(expression->subexpressions[counter], "");
      printf(separator);
    }
    printf(")");
  }
}


logical_expression * read_expression(char * input_string)
{
  long counter = 0;
  return read_expression(input_string, counter);
}


logical_expression * read_expression(char * input_string, long & counter)
{
  logical_expression * result = new logical_expression();
  long length =  string_length(input_string);
  while(1)
  {
    if (counter >= length)
    {
      break;
    }
    if (isspace(input_string[counter]))  
    {
      counter++;
      continue;
    }
    else if (input_string[counter] == '(')
    {
      counter++;
      read_word(input_string, counter, result->connective);
      read_subexpressions(input_string, counter, result->subexpressions);
      break;
    }
    else
    {
      read_word(input_string, counter, result->symbol);
      break;
    }
  }

  return result;
}


long read_subexpressions(char * input_string, long & counter, 
                         vector <logical_expression*> & subexpressions)
{
  long length =  string_length(input_string);
  while(1)
  {
    if (counter >= length)
    {
      printf("\nunexpected end of input\n");
      return 0;
    }
    if (isspace(input_string[counter]))    
    {
      counter++;
      continue;
    }
    if (input_string[counter] == ')') 
    {
      counter++;
      return 1;
    }
    else
    {
      logical_expression * expression = read_expression(input_string, counter);
      subexpressions.push_back(expression);
    }
  }
}


void read_word(char * input_string, long & counter, char * target)
{
  unsigned long second_counter = 0;
  while (1)
  {
    if (counter >= (long) string_length(input_string))
    {
      break;
    }
    if ((isalpha(input_string[counter])) || (input_string[counter] == '_') ||
        (isdigit(input_string[counter])))
    {
      target[second_counter] = input_string[counter];
      counter++;
      second_counter++;
    }
    else if ((input_string[counter] == ')') || (isspace(input_string[counter])))
    {
      break;
    }
    else
    {
      printf("unexpected character %c\n", input_string[counter]);
      exit_function(0);
    } 
  }

  target[second_counter] = 0;
}


long valid_expression(logical_expression *  expression)
{
  if (strcmp(expression->symbol, "") != 0)
  {
    return valid_symbol(expression->symbol);
  }

  if ((strcasecmp(expression->connective, "if") == 0) ||
      (strcasecmp(expression->connective, "iff") == 0))
  {
    if (expression->subexpressions.size() != 2)
    {
      printf("error: connective \"%s\" with %li arguments\n", 
             expression->connective, expression->subexpressions.size());
      return 0;
      return 0;
    }
  }
  else   if (strcasecmp(expression->connective, "not") == 0)
  {
    if (expression->subexpressions.size() != 1)
    {
      printf("error: connective \"%s\" with %li arguments\n", 
             expression->connective, expression->subexpressions.size());
      return 0;
    }
  }
  else if ((strcasecmp(expression->connective, "and") != 0) &&
           (strcasecmp(expression->connective, "or") != 0) &&
           (strcasecmp(expression->connective, "xor") != 0))
  {
    printf("error: unknown connective %s\n", expression->connective);
    return 0;
  }

  unsigned long counter;
  for (counter = 0; counter < expression->subexpressions.size(); counter++)
  {
    if (valid_expression(expression->subexpressions[counter]) == 0)
    {
      return 0;
    }
  }

  return 1;
}


long valid_symbol(char * symbol)
{
  if ((symbol == 0) || (strlen(symbol) == 0))
  {
    return 0;
  }
  unsigned long counter;
  for (counter = 0; counter < strlen(symbol);  counter++)
  {
    if ((symbol[counter] != '_') && (! (isalpha(symbol[counter]))) && (! (isdigit(symbol[counter]))))
    {
      return 0;
    }       
  }

  return 1;
}


int exit_function(int value)
{


  exit(value);
  return value;
}


void check_true_false(logical_expression * knowledge_base, 
                      logical_expression * statement)
{
 
 // CODE
  symbol_list = extract_symbol(knowledge_base);

  vector<bool> mod,known;
  mod.resize(symbol_list.size(), false);
  known.resize(symbol_list.size(), false);
  known = getKnown(knowledge_base);
  mod = getModel(knowledge_base);


  int output = entailment_function(knowledge_base,statement,symbol_list,mod, known);
  if(output == 0){
        printf("\ndefinitely true.\n");

  }
  else if(output == 1){
    printf("\ndefinitely false.\n");

  }
  else if(output == 2){
    printf("\n possibly true, possibly false.\n");

  }
  else{
    printf("\nBoth true and false.\n");

  }

  FILE * fp = fopen("result.txt", "wb");
  if (fp == 0)
  {
    printf("something is wrong, cannot open result.txt for writing\n");
  }
  else
  {
    if(output == 0){
        fprintf(fp,"\ndefinitely true.\n");

  }
  else if(output == 1){
    fprintf(fp,"\ndefinitely false.\n");

  }
  else if(output == 2){
    fprintf(fp,"\n possibly true, possibly false.\n");

  }
  else{
    fprintf(fp,"\nBoth true and false.\n");

  }
    fclose(fp);
  }
}

//------------------------------Entailment function--------------------------------

int entailment_function(logical_expression * kb, logical_expression * statement,vector<string>symbol_vec, vector<bool> model, vector<bool> known){

  logical_expression * neg_statement = new logical_expression();
  strcpy(neg_statement->connective, "not");
  neg_statement->subexpressions.push_back(statement);

  vector<string> unknown_symbols;
  for(int k = 0; k < known.size(); k++){
    if(!known[k]){
      unknown_symbols.push_back(symbol_vec[k]);
    }
  }
  bool A, B;
  A = isentail(kb,statement,symbol_vec,model,known, unknown_symbols);
  B = isentail(kb,neg_statement,symbol_vec,model,known,unknown_symbols);

  if(A==true && B==false){return 0;} // confirm true 
  else if(A==false && B==true){return 1;} // confirm false
  else if(A==true && B==true){return 2;} // Possibly either true or false
  else{return 3;} // Both true and false
}


bool isentail(logical_expression * kb, logical_expression * query,vector<string> symbol_vec, vector<bool> model, vector<bool> known, vector<string> unknown_symbols){
  
  if(unknown_symbols.empty()){
    if(evaluate_phrase(kb,model)){return evaluate_phrase(query,model);}
    else{return true;}
  }
  else{

    string P = unknown_symbols.front();
    unknown_symbols.erase(unknown_symbols.begin());
    vector<string>::iterator it_change = find(symbol_vec.begin(),symbol_vec.end(),P);

    vector<string>::iterator it = find(symbol_vec.begin(),symbol_vec.end(),unknown_symbols.front());
    
    vector<bool> temptruemod;
    vector<bool> tempfalsemod;
    for(int k = 0; k <distance(symbol_vec.begin(),it_change); k++){
        temptruemod.push_back(model[k]);
        tempfalsemod.push_back(model[k]);
    }
    temptruemod.push_back(true);
    tempfalsemod.push_back(false);
    
    return(isentail(kb,query,symbol_vec,temptruemod,known,unknown_symbols) && isentail(kb,query,symbol_vec,tempfalsemod,known,unknown_symbols));
  }
}

bool evaluate_phrase(logical_expression * statement,vector<bool> model){

  if(!strcmp(statement->connective,"or")){

          int temp = 0;
          for(int sub_it = 0; sub_it < statement->subexpressions.size();sub_it++){
            temp = temp + (int)evaluate_phrase(statement->subexpressions[sub_it],model);
          }
          if(temp == 0){return false;}
          else{return true;}
  }
  else if(!strcmp(statement->connective, "and")){

          int temp = 0;
          for(int sub_it = 0; sub_it < statement->subexpressions.size();sub_it++){
            temp = temp + (int)evaluate_phrase(statement->subexpressions[sub_it],model);
          }
          if(temp == statement->subexpressions.size()){return true;}
          else{return false;}
  }
  else if(!strcmp(statement->connective,"xor")){

          int temp = 0;
          for(int sub_it = 0; sub_it < statement->subexpressions.size();sub_it++){
            temp = temp + (int)evaluate_phrase(statement->subexpressions[sub_it],model);
          }
          if(temp == 1){return true;}
          else{return false;}
  }
  else if(!strcmp(statement->connective,"not")){
    
    if(evaluate_phrase(statement->subexpressions[0],model)==true){return false;}
    else{return true;}
  }
  else if(!strcmp(statement->connective,"if")){
        

          if(evaluate_phrase(statement->subexpressions[0],model) == true && evaluate_phrase(statement->subexpressions[1],model) == false){return false;}
          else{return true;}
  }
  else if(!strcmp(statement->connective,"iff")){

          

          if(evaluate_phrase(statement->subexpressions[0],model) == evaluate_phrase(statement->subexpressions[1],model)){return true;}
          else{return false;}
  }
  else{
          

          return(evaluate_symbol(statement->symbol,model));
  }
}


bool evaluate_symbol(char * symbol,vector<bool> model){
  vector<string>::iterator it = find(symbol_list.begin(),symbol_list.end(),symbol);
  return model[distance(symbol_list.begin(),it)];
}

vector<string> extract_symbol(logical_expression * statement){

  if(statement->subexpressions.size() == 0){
    vector<string> temp;
    temp.push_back(statement->symbol);
    return temp;
  }
  else{
    vector<string> tempvec;
    for(int sub_it = 0; sub_it < statement->subexpressions.size(); sub_it++){
      vector<string> all_temp;
      all_temp = extract_symbol(statement->subexpressions[sub_it]);

      for(int all_it = 0; all_it < all_temp.size();all_it++){
        tempvec.push_back(all_temp[all_it]);
      }
      all_temp.clear();
    }
    sort(tempvec.begin(),tempvec.end());
    tempvec.erase(unique(tempvec.begin(),tempvec.end() ),tempvec.end() );
    return tempvec;
  }
}

vector<bool> getKnown(logical_expression * statement){
  vector<bool> known;
  known.resize(symbol_list.size(), false);

  for(int k = 0;k< statement->subexpressions.size();k++){
    if(strcmp(statement->subexpressions[k]->symbol,"") !=0){
      vector<string>::iterator it = find(symbol_list.begin(),symbol_list.end(),statement->subexpressions[k]->symbol);
      
      known[distance(symbol_list.begin(),it)] = true;
    }
    else if(strcmp(statement->subexpressions[k]->connective,"not") == 0){
      if(strcmp(statement->subexpressions[k]->subexpressions[0]->symbol,"") != 0){
        vector<string>::iterator it = find(symbol_list.begin(),symbol_list.end(),statement->subexpressions[k]->subexpressions[0]->symbol);
        
        known[distance(symbol_list.begin(),it)] = true;
      }
    }
  }
  return known;
}

vector<bool> getModel(logical_expression * statement){
  vector<bool> mod;
  mod.resize(symbol_list.size(), false);

  for(int k = 0;k< statement->subexpressions.size();k++){
    if(strcmp(statement->subexpressions[k]->symbol,"")){
      

      vector<string>::iterator it = find(symbol_list.begin(),symbol_list.end(),statement->subexpressions[k]->symbol);
      mod[distance(symbol_list.begin(),it)] = true;
      
    }
    else if(strcmp(statement->subexpressions[k]->connective,"not")){
      if(!strcmp(statement->subexpressions[k]->subexpressions[0]->symbol,"")){
        vector<string>::iterator it = find(symbol_list.begin(),symbol_list.end(),statement->subexpressions[k]->subexpressions[0]->symbol);
        
        mod[distance(symbol_list.begin(),it)] = false;
        
      }
    }
  }
  return mod;
}




int main(int argc, char ** argv)
{
  char ** command_line = argv;
  if(argc != 4) 
  {
 
 // Takes 2 argument

    printf("Usage: %s [wumpus-rules-file] [additional-knowledge-file] [input_file]\n", command_line [0]);
    return exit_function(0);
  }

  char buffer[200];
  char * input;
  FILE * input_file;

 // READING WUMPUS RULE

  input = command_line[1];
  input_file = fopen(input, "rb");
  if (input_file == 0)
  {
    printf("failed to open file %s\n", input);
    return exit_function(0);
  }

  printf("Loading wumpus rules...\n");
  logical_expression * knowledge_base = new logical_expression();
  strcpy(knowledge_base->connective, "and");
  while(fgets(buffer, 200, input_file) != NULL)
  {
  
    if ((buffer[0] == '#') || (buffer[0] == 0) || (buffer[0] == 13) || (buffer[0] == 10))
    {
      continue;
    }
    logical_expression * subexpression = read_expression(buffer);
    knowledge_base->subexpressions.push_back(subexpression);
  }
  fclose(input_file);
  

// READING  ADDITIONAL KNOWLEDGE INFORMATION

  input = command_line[2];
  input_file = fopen(input, "rb");
  if (input_file == 0)
  {
    printf("failed to open file %s\n", input);
    return exit_function(0);
  }

  printf("Loading additional knowledge...\n");
  while(fgets(buffer, 200, input_file) != NULL)
  {
   
    if ((buffer[0] == '#') || (buffer[0] == 0) || (buffer[0] == 13) || (buffer[0] == 10))
    {
      continue;
    }
    logical_expression * subexpression = read_expression(buffer);
    knowledge_base->subexpressions.push_back(subexpression);
  }
  fclose(input_file);
  
  if (valid_expression(knowledge_base) == 0)
  {
    printf("invalid knowledge base\n");
    return exit_function(0);
  }

  print_expression(knowledge_base, "\n");

  // READING STATEMENT WHOSE ENTAILMENT WE WANT TO DETERMINE

  input = command_line[3];
  input_file = fopen(input, "rb");
  if (input_file == 0)
  {
    printf("failed to open file %s\n", input);
    return exit_function(0);
  }

  printf("\n\nLoading statement...\n");
  fgets(buffer, 200, input_file);
  fclose(input_file);

  logical_expression * statement = read_expression(buffer);
  if (valid_expression(statement) == 0)
  {
    printf("invalid statement\n");
    return exit_function(0);
  }

  print_expression(statement, "");

  check_true_false(knowledge_base, statement);

  delete knowledge_base;
  delete statement;
  exit_function(1);
}
