/*
 * Created for build automation.
 */
package org.reactome.fi;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Headless CLI entry point for the FI network build. The pipeline steps
 * (FINetworkBuilder, ReactomeDatabaseModifier, ...) are plain classes with
 * parameterless methods that were previously run one at a time as JUnit @Test
 * methods from Eclipse. This class reflectively invokes them in sequence so the
 * whole pipeline (or a single step) can be run from one command, with combined
 * logging and a hard stop on the first failure.
 *
 * Usage:
 *   mvn compile exec:java -Dexec.mainClass=org.reactome.fi.PipelineRunner
 *   mvn compile exec:java -Dexec.mainClass=org.reactome.fi.PipelineRunner \
 *       -Dexec.args="org.reactome.fi.FINetworkBuilder.buildFIDb"
 *
 * The full 6-step FINetworkBuilder sequence has two manual checkpoints baked into
 * it (see the class javadoc on FINetworkBuilder: a human should eyeball the
 * converted pathway projects before dumpPathwayDBs, and buildFIDb depends on
 * RF prediction files that only exist after a manual training run + threshold
 * pick). Running all 6 steps back-to-back with no args used to skip both -
 * scripts/run_pipeline.sh --stage1/--stage2/--stage3 is the intended entry point
 * since it stops and prints instructions at each checkpoint; this class's no-arg
 * default is deliberately just the first, RF-independent stage so that directly
 * invoking PipelineRunner without the wrapper can't skip the checkpoints either.
 */
public class PipelineRunner {

    private static final List<String> DEFAULT_STEPS = Arrays.asList(
            "org.reactome.fi.FINetworkBuilder.prepareMappingFiles",
            "org.reactome.fi.FINetworkBuilder.convertPathwayDBs");

    public static void main(String[] args) throws Exception {
        installLogTee();

        List<String> steps = args.length == 0 ? DEFAULT_STEPS : Arrays.asList(args);
        Map<String, Object> instanceCache = new HashMap<String, Object>();

        for (String step : steps) {
            int lastDot = step.lastIndexOf('.');
            if (lastDot < 0) {
                System.err.println("STEP FAILED: " + step
                        + " (expected fully.qualified.ClassName.methodName)");
                System.exit(1);
            }
            String className = step.substring(0, lastDot);
            String methodName = step.substring(lastDot + 1);
            long start = System.currentTimeMillis();
            try {
                Object instance = instanceCache.get(className);
                if (instance == null) {
                    Class<?> clazz = Class.forName(className);
                    instance = clazz.getDeclaredConstructor().newInstance();
                    instanceCache.put(className, instance);
                }
                Method method = instance.getClass().getMethod(methodName);
                method.invoke(instance);
                double seconds = (System.currentTimeMillis() - start) / 1000.0;
                System.out.println(String.format("STEP OK: %s (%.1fs)", step, seconds));
            }
            catch (InvocationTargetException e) {
                System.err.println("STEP FAILED: " + step);
                e.getTargetException().printStackTrace();
                System.exit(1);
            }
            catch (Exception e) {
                System.err.println("STEP FAILED: " + step);
                e.printStackTrace();
                System.exit(1);
            }
        }
        System.out.println("All steps completed successfully.");
    }

    /**
     * Tees System.out/System.err into a combined per-run log file under logs/,
     * installed before any pipeline class is instantiated so log4j's
     * ConsoleAppender (configured by each class's own constructor via
     * PropertyConfigurator.configure) binds to the tee'd stream.
     */
    private static void installLogTee() throws IOException {
        File logsDir = new File("logs");
        if (!logsDir.exists())
            logsDir.mkdirs();
        String timestamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date());
        File logFile = new File(logsDir, "build_" + timestamp + ".log");
        FileOutputStream fos = new FileOutputStream(logFile, true);
        System.out.println("Combined log: " + logFile.getPath());
        System.setOut(new PrintStream(new TeeOutputStream(System.out, fos), true));
        System.setErr(new PrintStream(new TeeOutputStream(System.err, fos), true));
    }

    private static class TeeOutputStream extends OutputStream {
        private final OutputStream first;
        private final OutputStream second;

        TeeOutputStream(OutputStream first, OutputStream second) {
            this.first = first;
            this.second = second;
        }

        @Override
        public void write(int b) throws IOException {
            first.write(b);
            second.write(b);
        }

        @Override
        public void write(byte[] b, int off, int len) throws IOException {
            first.write(b, off, len);
            second.write(b, off, len);
        }

        @Override
        public void flush() throws IOException {
            first.flush();
            second.flush();
        }
    }
}
