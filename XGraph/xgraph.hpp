#ifndef XGRAPH_HPP
#define XGRAPH_HPP

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

template <typename T>
class Element {
public:
	Element() = default;

	Element(const std::string& label) : label_(label) {};

	Element(const std::string& label, const T& content) : label_(label), content_(content) {};

	inline std::string label() const;

	inline void relabel(const std::string& newLabel);

	inline T content() const;

	inline void set_content(const T& newContent);

	inline std::unordered_map<std::string, Element<T>*> contactors() const;

	inline void insert(const Element& another);

	inline void remove(const std::string label);

	inline void operator=(const Element& another);

	inline bool operator==(const Element& another) const;

	inline int degree() const;

private:
	std::string label_;
	std::unordered_map<std::string, Element<T>*> contactors_;
	T content_;
};

#endif
