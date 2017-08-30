var channel = Backbone.Radio.channel('training');

//
// Generic Models
//

IterableCollection = Backbone.Collection.extend({
  defaults: {
    "completed": false,
    "selected": false,
  },
  get_active: function() {
    if(!this.findWhere({"selected": true})) {
      var m = this.at(0);
      m.set('selected', true);
    }
    return this.findWhere({"selected": true});
  },
  get_next: function() {
    var m = this.findWhere({"selected": true});
    if(!m) { return false; }

    //-- assume if you're proceeding, you've completed the past unit
    m.set("completed", true);

    var next = this.at(this.indexOf(m)+1);
    if(!next) { return false; }

    this.invoke('set', {"selected": false})
    next.set("selected", true);
    return next;
  }
})

TrainingInstruction = Backbone.RelationalModel.extend({
  /* A line of text we want to call attention to */
  defaults: {
    "text": "",
    "completed": false
  }
});

TrainingInstructionCollection = IterableCollection.extend({
  /* A Step may have sub-sections to progress through (eg. text instructions) */
  model: TrainingInstruction,
});

TrainingStep = Backbone.RelationalModel.extend({
  /* A specific Instruction or UI challenge within a Module */
  defaults: {
    "order": 0,
    "name": "",
    "description": "",
    "completed": false,
    "selected": false
  },
  relations: [{
    type: 'HasMany',
    key: 'instructions',
    relatedModel: TrainingInstruction,
    collectionType: TrainingInstructionCollection,
  }],
});

TrainingStepCollection = IterableCollection.extend({
  /* The many steps that compose a Module */
  model: TrainingStep,
});

TrainingModule = Backbone.RelationalModel.extend({
  /* A section (Level) within training that relates to a specific goal */
  defaults: {
    "name": "",
    "selected": false,
    "completed": false,
    "level": null  // Used to associated with saving progress and general order
  },
  relations: [{
    type: 'HasMany',
    key: 'steps',
    relatedModel: TrainingStep,
    collectionType: TrainingStepCollection,
  }],
});

TrainingModuleCollection = IterableCollection.extend({
  /* The different sections (Levels) within training */
  model: TrainingModule,
  url: function() {
    return '/api/training/'+this.task_type+'/'
  },
  initialize: function(models, options) {
    this.task_type = options.task_type;
  }
});

TrainingLevel = Backbone.RelationalModel.extend({
  /* Directly associated with a Module, but is the database way
   * of storing progress */
  defaults: {
    "level": null,
    "last_created": null,
    "completions": 0
  }
});

TrainingLevelCollection = IterableCollection.extend({
  /* All of user's past training completions */
  model: TrainingLevel
});

TrainingTask = Backbone.RelationalModel.extend({
  /* Representative of which Task the training will prepare the user for */
  defaults: {
    "selected": false,
    "task": ""
  },
  relations: [{
    type: 'HasMany',
    key: 'progress',
    relatedModel: TrainingLevel,
    collectionType: TrainingLevelCollection,
  }]
}),

TrainingTaskCollection = IterableCollection.extend({
  /* List of all available training types there are */
  url: '/api/training/',
  model: TrainingTask,
});

//
// Generic Views
//

// -- Progress Views
TrainingStepProgressItem = Backbone.Marionette.View.extend({
  template: '#training-step-progress-item-template',
  tagName: 'li',
  className: 'list-inline-item',
  onRender: function() {
    if(this.model.get('selected')){
      this.$el.css({'color': 'blue'});
    }
  },
  modelEvents: {
    'change:selected': 'render',
    'change:completed': 'render'
  }
});

TrainingStepProgress = Backbone.Marionette.CollectionView.extend({
  childView: TrainingStepProgressItem,
  tagName: 'ul',
  className: 'list-inline text-center'
});

TrainingModuleProgressItem = Backbone.Marionette.View.extend({
  template: _.template('<%= name %>'),
  tagName: 'li',
  className: function() {
    return 'breadcrumb-item'+ (this.model.get('selected') ? ' active' : '');
  }
});

TrainingModuleProgress = Backbone.Marionette.CollectionView.extend({
  childView: TrainingModuleProgressItem,
  tagName: 'ol',
  className: 'breadcrumb'
});

TrainingModuleNavigation = Backbone.Marionette.View.extend({
  template: '#training-navigation-template',

  regions: {
    'steps': '#progress-steps'
  },

  onRender: function() {
    this.showChildView('steps', new TrainingModuleProgress({'collection': this.collection}));
  },

  collectionEvents: {
    'change:selected': 'render',
  },
});

//--
TrainingFooterView = Backbone.Marionette.View.extend({
  /* this.model = TrainingModule
  * this.collection = None
  */
  template: '#training-footer-template',
  templateContext: function() {
    var counts = this.model.get("instructions").pluck("completed");
    return {'progress': (_.compact(counts).length / counts.length)*100}
  },
  className: 'row justify-content-center my-3',

  initialize: function() {
    this.listenTo(this.model.get('instructions'), 'change:completed', this.render);
  },

  events: {
    'mousedown a': function() { channel.trigger('training:next:instruction'); }
  },
})

//-- Inheritable parents views
TrainingStepInstructionView = Backbone.Marionette.View.extend({
  /* Display for a TrainingInstruction instance */
  template: _.template('<%= text %>')
})

TrainingStepTextView = Backbone.Marionette.View.extend({
  /* this.model = TrainingStep
   * this.collection = None
   */
  template: '#training-text-template',
  className: 'row',

  regions: {
    'instructions': '#instructions'
  },

  initialize: function() {
    var self = this;
    this.listenTo(channel, 'training:next:instruction', function() {
      var m = self.model.get('instructions').get_next();
      if(m) {
        //-- If there are more Instructions to go through
        self.showChildView('instructions', new TrainingStepInstructionView({'model': m}));
      } else {
        //-- Bubble up to go to next Step
        channel.trigger('training:next:step');
      }
    });
  },

  onRender: function() {
    var instruction = this.model.get('instructions').get_active();
    this.showChildView('instructions', new TrainingStepInstructionView({'model': instruction}));
  }

})

TrainingStepView = Backbone.Marionette.View.extend({
  /* this.model = TrainingStep
   * this.collection = None
   */
  template: '#training-module-step-template',
  className: 'row',

  regions: {
    'text': '#training-text',
    'action': '#training-action',
    'footer': '#training-footer'
  },

  initialize: function() {
    var self = this;
    this.listenTo(channel, 'training:next:step', function() {
      var m = self.model.collection.get_next();
      if(m) {
        //-- If there are more Steps to go through
        //self.render()
      } else {
        //-- Bubble up to go to next Module
        channel.trigger('training:next:module');
      }
    });
  },

});

TrainingModuleView = Backbone.Marionette.View.extend({
  /* this.model = TrainingModule
  * this.collection = None
  */
  template: '#training-module-template',
  className: 'row',

  regions: {
    'progress': '#module-progress',
    'step': '#step-action'
  },

  initialize: function() {
    var self = this;
    this.listenTo(channel, 'training:next:module', function() {
      console.log('Next Module', self);
      var m = self.model.collection.get_next();
      if(m) {
        //-- Move onto the next Module
        self.model = m;
        self.render();
      }
    });
  },
});


TrainingTaskView = Backbone.Marionette.View.extend({
  /*
   * this.model = None
   * this.collection = TrainingModuleList
   */
  template: '#training-task-template',

  className: 'row',

  regions: {
    'navigation': '#training-navigation',
    'message': '#training-message-center',
    'module': '#training-module',
  },

  collectionEvents: {
    'sync': 'render',
  },

})

TrainingView = Backbone.Marionette.View.extend({
  /* Dedicated page for viewing training info and selecting additional training
   * this.model = None
   * this.collection = TrainingTaskCollection
   */
  template: '#training-template',

  regions: {
    'selection': '#selection',
    'content': '#training'
  },

  initialize: function() {
    this.collection = new TrainingTaskCollection();
    this.collection.fetch();
  },

  collectionEvents: {
    sync: 'render'
  },

  onRender: function() {
    if(this.collection.length) {
      var task = this.collection.get_active();

      if(task.get("task") == "r") {
        this.showChildView('content', new RETrainingTaskView());
      }

      // if(task.get("task") == "e") {
      //   this.showChildView('content', new NERTrainingTaskView());
      // }

    }
  }

})
