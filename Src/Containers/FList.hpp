#ifndef FLIST_HPP
#define FLIST_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"


/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FList
 * Please read the license
 *
 * This class is a linked list container.
 * It is a very basic list to enable strong performance.
 *
 * Please refere to unit test flistUTest.cpp to know more.
 */
template< class Object >
class FList {
        /** A list node */
        struct Node {
                Object target;	//< Object of the node
                Node* next;	//< Next node
        };

        Node* root; //< Root node, NULL if size is 0
        int size;		 //< Elements in the list

        /**
        * Copy a list into current object
        * The current list has to be empty when this function is called
        */
        void copy(const FList& other){
                const Node* FRestrict  otherRoot = other.root;
                Node * FRestrict * myRoot = &this->root;
                while(otherRoot){
                        (*myRoot) = new Node;
                        (*myRoot)->target = otherRoot->target;

                        myRoot = &(*myRoot)->next;
                        otherRoot = otherRoot->next;
                }
                *myRoot = 0;
                this->size = other.size;
        }

public:
        /** Constructor (of an empty list) */
        FList() : root(0) , size(0) {
        }

        /** Desctructor */
        virtual ~FList(){
                clear();
        }

        /**
        * Copy operator
        * This will clear the current list before copying
        * @param other the source list
        * @return the current list as a reference
        */
        FList& operator=(const FList& other){
            if(this != &other){
                clear();
                copy(other);
            }
            return *this;
        }

        /**
        * Copy constructor
        * @param other the source/original list
        */
        FList(const FList& other): root(0) , size(0)  {
                copy(other);
        }

        /**
        * To clear the list
        * Size is 0 after calling this function
        */
        void clear(){
                while(this->root){
                        Node*const FRestrict next = this->root->next;
                        delete this->root;
                        this->root = next;
                }
                this->size = 0;
        }

        /**
        * Push an element in the head of the list
        * @param inObject the object to insert
        */
        void push(const Object& inObject){
                Node* newNode   = new Node;
                newNode->target = inObject;
                newNode->next 	= this->root;

                this->root 	= newNode;
                ++this->size;
        }


        /**
        * To get head value (last pushed value)
        * if size == 0 then defaultValue is returned
        * the list is unchanged after this function
        * @param defaultValue as the returned value in case size == 0, equal Object() if no param as passed
        * @return first value if exists or defaultValue otherwise
        */
        Object& head(Object& defaultValue = Object()){
            if(this->size) return this->root->target;
            else return defaultValue;
        }

        /**
        * To get head value as const
        * if size == 0 then defaultValue is return
        * the list is unchanged after this function
        * @param defaultValue as the returned value in case size == 0, equal Object() if no param as passed
        * @return first value if exists or defaultValue otherwise
        */
        const Object& head(const Object& defaultValue = Object()) const {
            if(this->size) return this->root->target;
            else return defaultValue;
        }

        /**
        * To get the head value and remove it from the list
        * @return first value
        * @warning you must check the list's size before calling this function!
        */
        Object pop(){
            --this->size;
            Node* newNode   = this->root;
            this->root      = this->root->next;

            Object value    = newNode->target;
            delete newNode;

            return value;
        }

        /**
        * To get the number of elements in the list
        * @return size
        */
        int getSize() const{
                return this->size;
        }

        /**
          * This iterator allow accessing list's elements
          * If you change the target list pointed by an iterator
          * you cannot used the iterator any more.
          * <code>
          * FList<int> values; <br>
          * // inserting ... <br>
          * FList<int>::BasicIterator iter(values); <br>
          * while( iter.hasNotFinished() ){ <br>
          *     iter.data() += 1; <br>
          *     iter.gotoNext(); <br>
          * } <br>
          * </code>
          */
        class BasicIterator {
        private:
            Node* iter; //< current node on the list

        public:
            /**
              * Constructor needs the target list
              * @param the list to iterate on
              */
            BasicIterator(FList& list) : iter(list.root){
            }

            /** To gotoNext on the list */
            void gotoNext(){
                if(this->iter){
                    this->iter = this->iter->next;
                    if(this->iter) Prefetch_Write(this->iter->next);
                }
            }

            /**
            * Current pointed value
            * current iterator must be valide (hasNotFinished()) to use this function
            */
            Object& data(){
                return this->iter->target;
            }

            /**
            * Current pointed value
            * current iterator must be valide (hasNotFinished()) to use this function
            */
            const Object& data() const{
                return this->iter->target;
            }

            /**
            * To know if an iterator is at the end of the list
            * @return true if the current iterator can gotoNext and access to value, else false
            */
            bool hasNotFinished() const{
                return iter;
            }

        };

        /**
          * This iterator allow accessing list's elements
          * If you change the target list pointed by an iterator
          * you cannot used the iterator any more.
          * <code>
          * FList<int> values; <br>
          * // inserting ... <br>
          * FList<int>::ConstBasicIterator iter(values); <br>
          * while( iter.hasNotFinished() ){ <br>
          *     iter.data() += 1; <br>
          *     iter.gotoNext(); <br>
          * } <br>
          * </code>
          */
        class ConstBasicIterator {
        private:
            const Node* iter; //< current node on the list

        public:
            /**
              * Constructor needs the target list
              * @param the list to iterate on
              */
            ConstBasicIterator(const FList& list) : iter(list.root){
            }

            /** to gotoNext on the list */
            void gotoNext(){
                if(this->iter){
                    this->iter = this->iter->next;
                    if(this->iter) Prefetch_Read(this->iter->next);
                }
            }

            /**
            * Current pointed value
            * current iterator must be valide (hasNotFinished()) to use this function
            */
            Object data(){
                return this->iter->target;
            }

            /**
            * Current pointed value
            * current iterator must be valide (hasNotFinished()) to use this function
            */
            const Object& data() const{
                return this->iter->target;
            }

            /**
            * To know if an iterator is at the end of the list
            * @return true if the current iterator can gotoNext and access to value, else false
            */
            bool hasNotFinished() const{
                return iter;
            }
        };

};

#endif //FLIST_HPP
// [--LICENSE--]
