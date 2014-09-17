"""
markov.py: Markov chain model library

Markov.py contains an implementation of Markov chain models of arbitrary order 
for sequences of characters. The main class of the library is the MarkovModel 
class, which represents an nth-order Markov chain model. It can be trained on 
strings, serialize to/from streams, and, once trained, compute the encoding 
costs of new strings.

Written for BME 205 Fall 2012
By Adam Novak
"""

import collections, itertools, math

class DefaultString(object):
    """
    Represents a string with a "default" character. The default character is
    returned when attempting to get characters before the beginning or after the
    end of the string that this object wraps.
    """
    def __init__(self, string, default):
        """
        Make a new DefaultString.
        string is the string that is being wrapped. string can actually be any 
        sequence type with 1-character strings as elements.
        default is a single-character string to return when a character outside 
        of the range of the backing string is requested.
        """
        
        # Holds the backing string
        self.string = string
        
        # Holds the default character
        self.default = default
        
    def __getitem__(self, index):
        """
        Return the corresponding character form the string if 0 <= index < 
        len(string), and default otherwise.
        index must be an integer or a slice object.
        If index is actually a slice object, this method computes the slice of 
        the virtual string and returns that, as a string.
        """
        
        if isinstance(index, slice):
            # Compute the slice
            # Pattern form StackOverflow
            # http://stackoverflow.com/questions/2936863/python-implementing-
            # slicing-in-getitem
            return "".join([self[i] for i in xrange(index.start or 0, 
                index.stop, index.step or 1)])
        else:
            # Must be an int
            if index < 0 or index >= len(self.string):
                return self.default
            else:
                return self.string[index]
                
    def __len__(self):
        """
        Return the length of the DefaultString.
        Really just the length of the backing string.
        """
        return len(self.string)

def count_kmers(string, k, start_stop=None, off_start=False, 
    end_behavior="none"):
    """
    Count the kmers of length k (int) in the string string.
    
    If start_stop is set to a single-character string, it is used as the
    start/stop character. If start_stop is not set, start and stop characters
    are not used and only kmers completely within the string are counted.
    
    if off_start is set, the first kmer counted is the one ending in the
    string's first character.
    
    end_behavior can be "none", "one", or "all". If set to "none", the last kmer
    counted is that ending in the last character of the string. This is the only
    permitted value if start_stop is not set. If end_behavior is "one", the last
    kmer counted is that ending in the first stop character. If end_behavior is
    "all", the last kmer counted is that starting with the string's last
    character.
    
    If k = 1, end_behavior "all" instead acts just like end_behavior "one",
    counting the kmer of just the stop character when it reaches the end of the
    string.
    
    Returns a collections.Counter object with all the kmer counts.
    
    If off_start is set, or end_behavior is not "none", start_stop must be set.
    
    The empty string is a special case. If start_stop is set, the empty string
    is counted as k start/stop characters. Otherwise, it is not counted.
    """
    
    # Check to make sure end_behavior is one of the available choices:
    if end_behavior not in ("none", "one", "all"):
        raise Exception("invalid end_behavior {}".format(end_behavior))
    
    
    # Check for attempting to go off the ends without a stat/stop character
    if (off_start or end_behavior != "none") and start_stop is None:
        raise Exception("Cannot go off of start or end without a start/stop "
            "character.")
    
    # Make end_behavior "all" work like end_behavior "one" for the special case
    # k=1. If order is 0 (and k is 1), and we're going past the end, we need to
    # count the last kmer as just the stop character. Otherwise a model with no
    # empty strings would never stop. This conflates the count of stop
    # characters and the count of empty strings for the 0-order model.
    if end_behavior == "all" and k == 1:
        end_behavior = "one"
    
    # Make a counter to hold all the kmer counts.
    counts = collections.Counter()
    
    if start_stop is not None:
        # Replace string with a DefaultString that fixes all the special cases
        # at the ends.
        string = DefaultString(string, start_stop)
    
    if len(string) == 0:
        # Handle the empty string as a special case: just a kmer of k start/stop
        # characters. If we sent it through the other branch, because the
        # start/stop characters are the same, we would not be able to
        # differentiate before-the-beginning and after-the-end kmers.
        if start_stop is not None:
            # If start_stop is set, count the empty string as k start/stop
            # characters.
            self.counts[start_stop * (self.order + 1)] += 1
        # If it's not set, don't count the empty string at all.
    else:
        # String is not empty.
        # Determine where to start and stop the kmers. If we want to stay within
        # the string, we start at 0 and stop at len(string) - k
        # kmer_start is where to start the first kmer
        # kmer_stop is where to *start* the *last* kmer
        if off_start:
            # Start the first kmer such that its last character is the first 
            # real character in the string.
            kmer_start = -(k - 1)
        else:
            # Start at the first character in the string
            kmer_start = 0
            
        if end_behavior == "none":
            # Start the last kmer so it ends right at the end of the string.
            kmer_stop = len(string) - k
        elif end_behavior == "one":
            # We want to count the last kmer as the one ending with the first
            # stop character.
            kmer_stop = len(string) - (k - 1)
        elif end_behavior == "all":
            # Count the last kmer as the one starting with the last real
            # character. We're (probably) using a higher-order model, so the
            # last kmer we count needs to include at least one real string
            # character. Otherwise it gets counted as an instance of the empty
            # string, which is not what we want in the higher-order case.
            kmer_stop = len(string) - 1
        
        # Now we know which kmers need counting. 
        
        # Count all the k-mers. Add 1 to the end since xrange wants a 1-past-
        # the-end and we just have the actual start position of the final kmer.
        for start in xrange(kmer_start, kmer_stop + 1):
            # This holds the k-mer we're counting as a string
            kmer = string[start:start + k]

            # Count the k-mer
            counts[kmer] += 1
            
        return counts

class MarkovModel(object):
    """
    Represents a Markov Chain model.
    """
    def __init__(self, order, alphabet, start_stop="="):
        """
        Make a new MarkovModel of the specified order. An Order k-1 Markov Chain 
        model operates on k-mers.
        
        order must be a non-negative integer.
        
        alphabet is a string of all possible characters in strings that the 
        model will be trained on or asked to compute encoding costs for.
        
        start_stop is an optional parameter specifying the default start/stop 
        character for the model. The start/stop character must not appear in any
        strings used to train the model, or in strings that the model is asked 
        to compute the encoding cost of.
        """
        # This holds the order of the model (int)
        self.order = order
        
        # This holds the alphabet as a set, for efficient x in alphabet 
        # determination
        self.alphabet = set(alphabet)
        
        # This holds counts of k-mers in a counter
        self.counts = collections.Counter()
        
        # This holds counts after pseudocount application
        # Filled in in add_pseudocount.
        # If we don't use pseudocouints, it just stays pointing at the Counter 
        # with the real counts. If we do use pseudocounts, it gets set to a new 
        # defaultdict of floats with the pseudocount-adjusted counts.
        self.adjusted_counts = self.counts
        
        # This holds a defaultdict of computed (log) probabilities for each
        # kmer. The default value will be -infinity, for "impossible" kmers.
        # Filled in in compute_probabilities.
        self.log_probabilities = None
        
        # This holds the start/stop character
        self.start_stop = start_stop
        
    def train(self, string):
        """
        Take all the k-mers in the string (including those starting before or 
        ending after the string) and use them to train the model.
        string can actually be any sequence type with 1-character strings as 
        elements. Characters not in the model's alphabet are stripped.
        
        The empty string is counted as a kmer of exactly k start/stop 
        characters. For a 0-order model, the k-mer starting 1 past the end is 
        counted, while for all higher orders, the last kmer counted is that 
        starting at the last character in the string.
        
        The empty string is handled appropriately: in the 0-order model, the 
        stop character will be counted once for every sequence, even the empty 
        ones; the empty sequence in the 0-order model has the same probability 
        as any sequence has of stopping. In higher-order models, the empty 
        string will be counted as the string of just start/stop characters, 
        which reverses appropriately and still represents the empty string 
        (though with a slightly different semantic interpretation). This kmer 
        will only be produced by the empty string; in higher-order models, the 
        probability of the empty string is distinct from the probability of 
        stopping a sequence while generating.
        
        The trained model can easily be reversed by reversing all the kmers and 
        keeping the counts. 
        """

        # Strip non-alphabet characters
        string = filter(lambda c: c in self.alphabet, string)
        
        # Add in counts of all kmers in the string
        self.counts += count_kmers(string, self.order + 1, 
            start_stop=self.start_stop, off_start=True, end_behavior="all")
    
    def set_pseudocount(self, pseudocount):
        """
        Add the specified constant pseudocount to each kmer possible in the 
        model's alphabet. Kmers starting before the beginning of the string,
        those ending after the end, and the kmer corresponding to the empty 
        string are computed appropriately based on the order of the model, and 
        also have the pseudocount added to them.
        
        pseudocount specifies the constant float pseudocount to add.
        
        If called more than once, only the pseudocount specified in the final 
        call will be honored.
        """        
        
        # Make a new defaultdict
        # The function passed here is called to intitalize each entry in the 
        # dict on demand. We have a function that, when called, produces a new
        # float with the value of count. Just float(count) won't work because 
        # you can't call that.
        self.adjusted_counts = collections.defaultdict(
            lambda: float(pseudocount))
            
        # Add in nonzero counts from the counter.
        for (kmer, count) in self.counts.iteritems():
            self.adjusted_counts[kmer] += count
        
        
    
    def __compute_probabilities_for_prefix(self, prefix):
        """
        Compute the probabilities for all kmers that start with the given 
        k-1-character prefix. The probability of a kmer is defined as its count 
        over the total count of all kmers that match its first k-1 characters. 
        Probabilities are stored in the internal "probabilities" dict, as log 
        probabilities.
        
        Special care must be taken when a prefix never appeared in the data, and
        no pseudocounts have been added. Probabilities for kmers where nothing
        in the prefix has a count are left at 0, or -infinity in log-probability
        terms.
        
        self.log_probabilities must not be None.
        """
        
        # This holds the sum of all counts for kmers with this prefix
        prefix_total = 0
        # Don't forget to count the start/stop character as a last letter
        for letter in itertools.chain(self.alphabet, self.start_stop):
            prefix_total += self.adjusted_counts[prefix + letter]
            
        if prefix_total == 0:
            # This whole prefix is unobserved, and not covered by pseudocounts
            # either. So what should we put for probabilities of kmers given
            # prefixes? Let's just leave them at their default value of
            # impossible.
            return
            
        # Now start dividing
        for letter in itertools.chain(self.alphabet, self.start_stop):
            # This holds the kmer we're finding the probability for, as a 
            # string
            kmer = prefix + letter

            if self.adjusted_counts[kmer] != 0:
                # This kmer actually exists!
                # logP(kmer|prefix) = kmer count / count of all things with  
                # that prefix
                self.log_probabilities[kmer] = math.log(
                    float(self.adjusted_counts[kmer]) / float(prefix_total), 2)
            # If a kmer has 0 count, it needs probability 0, or log probability
            # -infinity, which is the default.
        
    def compute_probabilities(self):
        """
        Compute the probabilities for each kmer. The probability of a kmer is 
        defined as its count over the total count of all kmers that match its 
        first k-1 characters. Probabilities are stored in the internal 
        "probabilities" dict, as log probabilities.
        
        Special care must be taken when a prefix never appeared in the data, 
        and no pseudocounts have been added. Probabilities for kmers where 
        nothing in the prefix has a count are set to 0.
        """
        
        # This is a defaultdict of floats of log-probabilities of kmers given
        # their prefixes. It starts out as None, so we fill it in here so
        # encoding_cost has something to work with. We make the default value
        # -infinity, so that "impossible" kmers are handled gracefully later (no
        # KeyErrors).
        self.log_probabilities = collections.defaultdict(lambda: float("-inf"))
        
        # The possible prefixes are tricky
        # First handle those of order alphabet chatracters
        # For order 0 we should just get the prefix ""
        # For order 0, this handles the empty string probability.
        for prefix in itertools.product(self.alphabet, repeat=self.order):
            self.__compute_probabilities_for_prefix("".join(prefix))
                
        # Now possible prefixes involving start/stop characters
        # This only does anything for order > 0, and handles the empty string 
        # probability there.
        # We only need to consider prefixes having start/stop characters in a 
        # contiguous block at the front. Nobody may ever ask for the probability
        # of anything coming after a stop.
        # The prefix can have up to k-1 = order start/stops, and must have
        # at least 1
        for num_start_stops in xrange(1, self.order + 1):
            # This holds a string with the specified number of start/stop 
            # characters.
            start_stops = self.start_stop * num_start_stops
            # For each combination of things from the alphabet that would 
            # bring us up to k-1 characters (possibly including "")...
            for combination in itertools.product(self.alphabet, 
                repeat=self.order - num_start_stops):
                
                # Use the concatenation as the prefix and find probabilities
                # for things in it
                self.__compute_probabilities_for_prefix(start_stops + 
                    "".join(combination))
                
    def encoding_cost(self, string):
        """
        Finds the encoding cost in bits for the given string. The encoding cost
        is the sum of -log base 2 of the probability of each character given the
        preceeding order characters.
        
        Assumes that compute_probabilities() has been called on the model.
        
        string is the string for which the encoding cost is to be computed. It 
        must not contain the start/stop characters. Any characters not in the 
        model's alphabet will be stripped from the string.
        
        Note that, if pseudocounts are not used, a string containing a kmer that
        does not appear in the model will have infinite encoding cost.
        
        Returns a tuple of the encoding cost and the number of characters in the
        string (that passed the alphabet filter). The end of the string counts 
        as a stop "character", so the empty string has one "character".
        """
    
        if self.log_probabilities is None:
            # They forgot to call calculate_probabilities()
            raise Exception("Cannot compute encoding cost without first "
                "calculating probabilities.")
    
        # Strip non-alphabet characters
        string = filter(lambda c: c in self.alphabet, string)
    
        # Count all the kmers in the string. Go off the start but only to the
        # first stop at the end, because kmers after the first stop character
        # mean nothing for probability. This variable gets a collections.Counter
        # of kmer counts in the string we're computing the encoding cost of.
        kmers = count_kmers(string, self.order + 1, start_stop=self.start_stop, 
            off_start=True, end_behavior="one")

        # This holds the running total encoding cost as a float
        cost = 0.0
        
        # This holds the number of kmers looked at. It ought to equal 
        # len(string) + 1, but we want to be safe, so we count them explicitly.
        kmer_count = 0
            
        # For each kmer kmer and its count count...
        for (kmer, count) in kmers.iteritems():
            # Subtract the kmer log probability times its count from the total,
            # since we're summing negative log probabilities. Don't need to
            # check for presence of key since this is a defaultdict.
            cost -= self.log_probabilities[kmer] * count
            
            # Count the number of times that kmer appeared.
            kmer_count += count
            
        return (cost, kmer_count)
        
    def write(self, stream):
        """
        Output all kmers the model has been trained on, in alphabetical order, 
        with their corresponding counts, one per line, to the given stream.
        
        If pseudocounts have been added, these will not appear in the output.
        """
        
        # This holds the item list that we need to sort
        item_list = self.counts.items()
        item_list.sort()
        
        # kmer gets each kmer string, and count gets its count
        for (kmer, count) in item_list:
            stream.write("{} {}\n".format(kmer, count))
    
    @classmethod        
    def read(cls, stream, alphabet, start_stop="="):
        """
        Reads a MarkovModel from the given stream and returns it. Auto-detects 
        order based on the length of the first kmer.
        
        This is a class method, so you call it directly from MarkovModel, not an
        instance.
        
        cls is the class on which this method is called. It will be MarkovModel 
        unless the method has been inherited by another class.
        
        stream is the stream to read from. All data from the stream before EOF 
        must constitute a valid model description, as written by 
        MarkovModel.write().
        
        alphabet is a string containing each acceptable character in kmers. If a
        kmer contains unacceptable characters, it will be skipped. alphabet is
        also the alphabet used to construct the MarkovModel object returned.
        
        start_stop is the character to be used as the start/stop character for 
        the model. This would be quite challenging to auto-detect from the model
        description. It must match the start/stop character used in the original
        model before serialization.
        """
        
        # This holds a set of all characters allowed in kmers.
        # Even though self.alphabet is a set, we need a new copy.
        acceptable_characters = set(alphabet)
        # The start/stop character is always acceptable.
        acceptable_characters.add(start_stop)
        
        # This holds the MarkovModel under construction.
        # The order passed here will be overwritten later.
        # We construct it using cls in case we've been inherited.
        model = cls(0, alphabet, start_stop)
        
        # This flag records if we've read the first line and set the order or 
        # not.
        order_set = False
        
        for line in stream:
            # kmer holds the kmer string, and count holds the count as a string
            (kmer, count) = line.split()
            
            # Skip kmers containing invalid characters
            # Compute set difference between kmer's characters and the 
            # acceptable characters.
            if len(set(kmer) - acceptable_characters) > 0:
                continue
            
            if not order_set:
                # Calculate the order
                model.order = len(kmer) - 1
                order_set = True
            
            # Add the count of the kmer to the model.    
            model.counts[kmer] += int(count)
        
        return model
            
