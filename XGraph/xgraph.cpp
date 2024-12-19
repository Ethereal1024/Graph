#include "xgraph.hpp"

template <typename T>
std::string Element<T>::label() const {
	return label_;
}

template <typename T>
void Element<T>::relabel(const std::string& newLabel) {
	label_ = newLabel;
}

template <typename T>
T Element<T>::content() const {
	return content_;
}

template <typename T>
void Element<T>::set_content(const T& newContent) {
	content_ = newContent;
}

template <typename T>
std::vector<Element<T>*> Element<T>::contactors() const {
	std::vector<Element<T>*> rst;
	for(const auto& pair : contactors_) {
		rst.push_back(pair.second);
	}
	return rst;
}

template <typename T>
void Element<T>::insert(const Element<T>& another) {
	contactors[another.label] = &another;
}

template <typename T>
void Element<T>::remove(const std::string label) {
	contactors_.erase(label);
}

template <typename T>
void Element<T>::operator=(const Element<T>& another) {
	label_ = another.label();
	contactors_ = another.contactors();
	content_ = another.content();
}

template <typename T>
bool Element<T>::operator==(const Element<T>& another) const {
	return label_ == another.label() && contactors_ == another.contactors() && content == another.content();
}

template <typename T>
int Element<T>::degree() const {
	return contactors_.size();
}


 
