package edu.ucsc.genome

/**
 * A collection of useful methods for setting up a command-line tool that hits
 * against the Sequence Graphs API. Mostly deals with the configuration of
 * logging, in order to quiet the built-in loggers so you can actually see your
 * tool's output.
 *
 * We can't just have the library unilaterally configure logging without the
 * application's permission; that's the sort of behavior on the part of *our*
 * dependencies that we're trying to work around.
 *
 * If you are using this class, you also probably need a log4j.properties file
 * that suppresses meaningless INFO messages from Spark and its many
 * dependencies scattered throughout namespace-world. Something like:
 *
 * {{{
 * log4j.rootLogger=WARN, stdout 
 * }}}
 */
object ConsoleUtil {
    /**
     * Make Parquet's built-in logger, which helpfully configures itself to log
     * all INFO messages to console in exactly the way that a good library
     * shouldn't, be quiet.
     */
    def quietParquet = {
        import java.util.logging._
        import parquet.Log

    
        // Parquet "helpfully" forcibly adds an INFO-level logger to console if
        // we don't configure its logger specifically (i.e. if we let its log
        // messages pass through to the root logger). It also overrides any
        // level we set on its logger, resetting it to INFO. See
        // <https://github.com/Parquet/parquet-mr/blob/master/parquet-
        // common/src/main/java/parquet/Log.java>
        
        // The solution is to add a warning-level console logger specifically
        // for Parquet, and tell it not to propagate messages up to the root
        // logger.
        
        // This could be done through the properties-file-based Java logging
        // configuration mechanism, but that would require telling Java how to
        // *find* the properties file, which in turn requires either setting a
        // JVM command-line option to a filesystem path or messing about with
        // the Java preferences API and/or resource streams. See
        // <http://stackoverflow.com/q/805701/402891>
        
        // Get the logger Parquet is going to use (before Parquet's static
        // logger initialization code can run)
        val parquetLogger = Logger.getLogger(
            classOf[parquet.Log].getPackage().getName())
        // Attach our own ConsoleHandler for it
        val consoleHandler = new ConsoleHandler()
        consoleHandler.setLevel(Level.WARNING)
        parquetLogger.addHandler(consoleHandler)
        // Tell it not to send messages up
        parquetLogger.setUseParentHandlers(false)
    
    }
    
    /**
     * Make log4j's root logger report only actually dangerous things, without
     * having to have a log4j.properties file packaged in your application's
     * jar. This is useful when you are lazy and don't want to add extra
     * resource files.
     */
    def quietLog4j = {
        import org.apache.log4j.{LogManager, Level}
        
        // Show only WARNs or worse. See
        // <http://stackoverflow.com/a/4598829/402891>
        LogManager.getRootLogger().setLevel(Level.WARN);
    }
}
