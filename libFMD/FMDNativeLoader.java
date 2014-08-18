package org.ga4gh;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * A class to load the FMD native library whenever the FMD jar is loaded.
 * Only works if it's packaged into a jar with "libfmd.so" in it, under the path
 * corresponding to the package in which this class is found.
 *
 * Based on:
 * Simple library class for working with JNI (Java Native Interface)
 * @see http://frommyplayground.com/how-to-load-native-jni-library-from-jar
 * @author Adam Heirnich <adam@adamh.cz>, http://www.adamh.cz
 * Licensed under the MIT license.
 */
public class FMDNativeLoader {
    static boolean loaded = false;
    static void load() {
        if(loaded) {
            // Don't try and load twice.
            return;
        }
        try {
        
            // Load the libraries which this library depends on. They must be
            // installed at the system level. See
            // <http://stackoverflow.com/q/5425034/402891>
            System.loadLibrary("boost_system");
            System.loadLibrary("boost_filesystem");
        
            // Read the library data out of our jar
        	InputStream libraryData = 
        	    FMDNativeLoader.class.getResourceAsStream("libfmd.so");
            if (libraryData == null) {
                // Complain we didn't find the library in our jar. See
                // <http://frommyplayground.com/how-to-load-native-jni-library-
                // from-jar/>
                throw new FileNotFoundException(
                    "Can't load libfmd.so from FMD jar!");
            }
            
            // Make a temporary directory to hold the library
            Path libraryDir = Files.createTempDirectory("libfmd");
            
            // Make a Path for the file we're loading
            Path library = libraryDir.resolve("libfmd.so");
            
            // Copy all the data over with the new NIO file copying method.
            Files.copy(libraryData, library);
            
            // Register the directory and then the library to be deleted.
            // deleteOnExit hooks happen in reverse order. See
            // <http://www.coderanch.com/t/278832/java-io/java/delete-directory-
            // VM-exits>
            libraryDir.toFile().deleteOnExit();
            library.toFile().deleteOnExit();
            
            // Load the library
        	System.load(library.toFile().getAbsolutePath());
        	
	    } catch(Exception e) {
	        // Catch any ordinary "checked" exceptions and rethrow them as 
	        // runtime exceptions, which are the only kind of exceptions static 
	        // blocks are allowed to throw. See 
            // <http://stackoverflow.com/a/15289277/402891>
	        throw new ExceptionInInitializerError(e);
	    }
	    
	    // Don't try to load again.
	    loaded = true;
        
    }
}
