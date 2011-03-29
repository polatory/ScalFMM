#ifndef FLIST_HPP
#define FLIST_HPP
// /!\ Please, you must read the license at the bottom of this page

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FList
 * Please read the license
 *
 * This class is a linked list container.
 * It is a very basic list to enable strong performance.
 *
 * Please refere to unit test flistUTest.cpp
 */
template< class Object >
class FList {
	/** A list node */
	struct Node {
		Object target;	//< Object of the node
		Node* next;	//< Next node
	};

	Node* root;		//< Root node, NULL if size is 0

        int size;		//< Elements in the list

        /**
        * Copy a list into current object
        * The current list has to be empty when this function is called
        */
	void copy(const FList& other){
		const Node* otherRoot = other.root;
		Node ** myRoot = &this->root;
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
		clear();
		copy(other);
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
			Node* next(this->root->next);
			delete this->root;
			this->root = next;
		}
		this->size = 0;
	}

	/**
	* Push an element in the front of the list
	*/
	void pushFront(Object& inObject){
		Node* newNode 	= new Node;
		newNode->target = inObject;
		newNode->next 	= this->root;

		this->root 	= newNode;
                ++this->size;
	}

        /**
        * Push an element in the front of the list
        * used when type is pointer and like [TYPE* const] object
        */
        void pushFront(Object const inObject){
                    Node* newNode           = new Node;
                    newNode->target         = inObject;
                    newNode->next           = this->root;

                    this->root  = newNode;
                    ++this->size;
        }

        /**
        * To get front value (last pushed value)
        * if size == 0 then defaultValue is returned
        * @param defaultValue as the returned value in case size == 0, equal Object() if no param as passed
        * @return first value if exists or defaultValue otherwise
        */
        Object& front(Object& defaultValue = Object()){
            if(this->size) return this->root->target;
            else return defaultValue;
        }

        /**
        * To get front value as const
        * if size == 0 then defaultValue is return
        * @param defaultValue as the returned value in case size == 0, equal Object() if no param as passed
        * @return first value if exists or defaultValue otherwise
        */
        const Object& front(const Object& defaultValue = Object()) const {
            if(this->size) return this->root->target;
            else return defaultValue;
        }

        /**
        * To get the front value and remove it from the list
        * @return first value
        * @warning you must check the list's size before calling this function!
        */
        Object popFront(){
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
          * while( iter.isValide() ){ <br>
          *     iter.value() += 1; <br>
          *     iter.progress(); <br>
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

            /** to progress on the list */
            void progress(){
                if(this->iter) this->iter = this->iter->next;
            }

            /**
            * Current pointed value
            * current iterator must be valide (isValide()) to use this function
            */
            Object& value(){
                return this->iter->target;
            }

            /**
            * Current pointed value
            * current iterator must be valide (isValide()) to use this function
            */
            const Object& value() const{
                return this->iter->target;
            }

            /**
            * To know if an iterator is at the end of the list
            * @return true if the current iterator can progress and access to value, else false
            */
            bool isValide() const{
                return iter;
            }

        };

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

            /** to progress on the list */
            void progress(){
                if(this->iter) this->iter = this->iter->next;
            }

            /**
            * Current pointed value
            * current iterator must be valide (isValide()) to use this function
            */
            Object value(){
                return this->iter->target;
            }

            /**
            * Current pointed value
            * current iterator must be valide (isValide()) to use this function
            */
            const Object& value() const{
                return this->iter->target;
            }

            /**
            * To know if an iterator is at the end of the list
            * @return true if the current iterator can progress and access to value, else false
            */
            bool isValide() const{
                return iter;
            }
        };

};

#endif //FLIST_HPP
// [--LICENSE--]
