"use client";

interface GeneErrorFallbackProps {
  id: string;
  error: string;
}

export default function GeneErrorFallback({ id, error }: GeneErrorFallbackProps) {
  return (
    <div className="container mx-auto py-20 px-4 text-center">
      <h2 className="text-2xl font-bold text-slate-300 mb-4">Dashboard Temporarily Unavailable</h2>
      <p className="text-slate-500 mb-8">
        The backend analytics service for Gene {id} is not responding. 
        <br />
        <span className="text-xs font-mono mt-2 block opacity-60">Error: {error}</span>
      </p>
      <button 
        onClick={() => window.location.reload()}
        className="px-6 py-2 bg-primary/10 border border-primary/20 text-primary rounded-lg hover:bg-primary/20 transition-all font-semibold"
      >
        Re-establish Connection
      </button>
    </div>
  );
}
