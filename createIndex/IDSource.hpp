/**
 * Represents a source for sequential IDs of the specified integer type.
 */
template <class T>
class IDSource {
public:
    /**
     * Create a new IDSource starting at the given ID.
     */
    IDSource(T start): nextID(start) {
        // Nothing to do
    }

    /**
     * Get the next available ID.
     */
    T next() {
        return nextID++;
    }

private:
    /**
     * Hold the next ID to return.
     */
    T nextID;
};
