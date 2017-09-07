var channel = Backbone.Radio.channel('training');

//
// Generic Models
//

IterableCollection = Backbone.Collection.extend({
  defaults: {
    "completed": false,
    "selected": false,
  },
  set_active: function(model) {
    this.invoke('set', {"selected": false})
    model.set("selected", true);
    return this.get_active();
  },
  get_active: function() {
    if(!this.findWhere({"selected": true})) {
      if(this.length >= 1) {
        var m = this.at(0);
        m.set('selected', true);
      } else {
        return false;
      }
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
    reverseRelation: {
      'key': 'module'
    }
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
  className: 'list-inline text-center mb-0'
});

TrainingModuleProgressItem = Backbone.Marionette.View.extend({
  template: _.template('<%= name %>'),
  tagName: 'li',
  className: function() {
    return 'breadcrumb-item'+ (this.model.get('selected') ? ' active' : '');
  },
  events: {
    'mousedown': function() {
      channel.trigger('training:goto:module', this.model);
    }
  }
});

TrainingModuleProgress = Backbone.Marionette.CollectionView.extend({
  childView: TrainingModuleProgressItem,
  tagName: 'ol',
  className: 'breadcrumb mb-0'
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
//-- Inheritable parents views
TrainingStepInstructionView = Backbone.Marionette.View.extend({
  /* Display for a TrainingInstruction instance
  * TrainingInstruction
  */
  template: '#training-instruction-text-template',
  templateContext: function() {
    return {'idx':  this.model.collection.indexOf(this.model)+1}
  },
  initialize: function() {
    if(this.model.get('training_data')) {
      channel.trigger('training:hide:action');
      channel.trigger('training:show:action', this.model.get('training_data'), this.model.get('training_rules'), this.model.collection.indexOf(this.model));
    } else {
      channel.trigger('training:hide:action');
    }
  }
})

TrainingStepTextView = Backbone.Marionette.View.extend({
  /* Text at the top half of the Insturction page, usually used to describe what to do,
   * and list the various instruction steps that need to be progressed through
   * this.model = TrainingStep
   * this.collection = None
   */
  template: '#training-text-template',
  className: 'row justify-content-center',

  ui: {
    'instructions': '#training-instructions'
  },

  regions: {
    'instructions': '#training-instructions'
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

    if(this.model.get('training_data')) {
      channel.trigger('training:hide:action');
      channel.trigger('training:show:action', this.model.get('training_data'), this.model.get('training_rules'))
    } else {
      channel.trigger('training:hide:action');
    }

  },

  onRender: function() {
    var instruction = this.model.get('instructions').get_active();
    if(instruction) {
      this.showChildView('instructions', new TrainingStepInstructionView({'model': instruction}));
    } else {
      this.ui.instructions.hide();
    }
  }

});

TrainingFooterView = Backbone.Marionette.View.extend({
  /* Displays the progression thorughout
  * this.model = TrainingModule
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
});


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

    this.listenTo(channel, 'training:show:action', function(training_data, training_rules, instruction_idx) {
      var step_idx = self.model.get('order');
      var module_idx = self.model.get('module').get('level');
      var instruction_idx = instruction_idx || null;
      this.showChildView('action', new RETrainingAction({'position': [module_idx, step_idx, instruction_idx], 'training_data': training_data, 'training_rules': training_rules}));
    });

    this.listenTo(channel, 'training:hide:action', function() {
      this.getRegion('action').empty();
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

    this.listenTo(channel, 'training:goto:module', function(module) {
      self.model = module.collection.set_active(module);
      self.render();
    });

  },
});


TrainingTaskView = Backbone.Marionette.View.extend({
  /* Base template for organizing everythin needed for going through a Task
   * this.model = None
   * this.collection = TrainingModuleList
   */
  template: '#training-task-template',

  className: 'row',

  regions: {
    'module-navigation': '#training-module-navigation',
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
   *
   * Training: Task (RE or NER) >> Module >> Step >> Instruction
   *
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
