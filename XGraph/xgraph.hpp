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

	inline std::vector<Element<T>*> contactors() const;

	inline void insert(const Element<T>& another);

	inline void remove(const std::string label);

	inline void operator=(const Element<T>& another);

	inline bool operator==(const Element<T>& another) const;

	inline int degree() const;

	virtual std::vector<Element<T>*> neighbors();

	virtual std::vector<std::string> neighbors_label();

private:
	std::string label_;
	std::unordered_map<std::string, Element<T>*> contactors_;
	T content_;
};

class Edge;

template <typename T>
class Vertex : private Element<T> {
public:
	void connect_to(const std::string& label);

	void disconnect_from(const std::string& label);

	std::vector<Vertex<T>*> neighbors() const override;

	std::vector<std::string> neighbors_label() const override;
private:
	using Element<T>::insert;
	
	using Element<T>::remove;
};

#endif
